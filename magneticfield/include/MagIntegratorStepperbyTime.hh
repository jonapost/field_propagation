// Nystrom stepper implemenation by Jason Suagee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 27 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license

#ifndef MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_
#define MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_


#include "G4Mag_EqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ThreeVector.hh"

#include "G4CachedMagneticField.hh"

#include <assert.h>
#include <vector>

using namespace std;

#include "isTracking.hh"

#define NO_STATE_VARIABLES 12

// Not all of these are needed (they are used in StepTracker.cc).
// Here we only use POSITION_SLOT:
// (Everything is offset by 2, one for time and one for arclength.)
#define POSITION_SLOT 2
#define MOMENTUM_SLOT 5
#define RHS_SLOT 8
#define BUFFER_COLUMN_LEN 11


template <class BaseStepper>

class MagIntegratorStepper_byTime : public BaseStepper {
public:
   inline
   MagIntegratorStepper_byTime(G4Mag_EqRhs *EquationRhs,
                                 G4int numberOfVariables = 6,
                                 G4int numStateVariables = 12);
   virtual ~MagIntegratorStepper_byTime();

   inline void ComputeRightHandSide(const G4double yInput[], G4double dydx[]);



   inline void Stepper(const G4double yInput[],
                       const G4double dydx[],
                             G4double hstep,
                             G4double yOutput[],
                             G4double yError [] );

   inline G4double DistChord() const;

private:
   // Needed because Stepper() takes yInput[] and dydx[] as const input
   // but we have to do some scaling modifications in Stepper():
   G4double yIn[NO_STATE_VARIABLES],
            dydx_copy[NO_STATE_VARIABLES];

   // Possible future use:
   //G4double cached_dydx[NO_STATE_VARIABLES],
   //         last_function_evaluation[NO_STATE_VARIABLES];

   G4double nextFunctionEvaluation[NO_STATE_VARIABLES];

   G4Mag_EqRhs  *m_fEq;
};


template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::ComputeRightHandSide(
                                             const G4double yInput[], G4double dydx[] ) {

   for (int i = 0; i < NO_STATE_VARIABLES; i ++)
      yIn[i] = yInput[i];
   for (int i = 3; i < 6; i ++)
      yIn[i] *= 1. / m_fEq -> FMass();

   BaseStepper::ComputeRightHandSide(yIn, dydx);

   for (int i = 3; i < 6; i ++) {      // may want to make this from i = 0, (but really
                                       // doesn't matter since we're only using this
                                       // for Nystrom Steppers.
      dydx[i] *= m_fEq -> FMass();
   }
   // Always feed this template class a Mag_UsualEqRhs_IntegrateByTime as EquationRhs.
}

template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::Stepper(const G4double yInput[],
      const G4double dydx[],
            G4double hstep,
            G4double yOutput[],
            G4double yError [] ) {

   // Have to copy because yInput and dydx are constant in the function signature...
   for (int i = 0; i < 6; i ++){
      yIn[i] = yInput[i];
      dydx_copy[i] = dydx[i];
   }
   // ...because now we have to convert to velocity coordinates:
   for (int i = 3; i < 6; i ++){
      yIn[i] *= 1. / m_fEq -> FMass();
      dydx_copy[i] *= 1. / m_fEq -> FMass();
   }

   /////// Getting velocity for purposes of rescaling the step length:
   G4double velocity = G4ThreeVector( yIn[3], yIn[4], yIn[5] ).mag();

   // Within stepper call, convert hstep (which is in units of arclength)
   // to hstep / velocity (which is units of time).
   BaseStepper::Stepper( yIn, dydx_copy, hstep / velocity, yOutput, yError );



#ifdef TRACKING
   if ( BaseStepper::mTracker -> get_within_AdvanceChordLimited() ) {
      // Within AdvancedChordLimited, so we want to record.

      if ( BaseStepper::mTracker -> isArmed() ) { // If armed then store step result.

         // We need the next RHS function evaluation (for the right endpoint of
         // step interval). We need this because it currently is not stored by
         // the stepper (FSAL?)

         // Also, because we are currently in velocity coordinates, we just call
         // the BaseStepper::ComputeRightHandSide() method.
         BaseStepper::ComputeRightHandSide(yOutput, nextFunctionEvaluation);

         // Now we have to undo the effect of Mag_UsualEqRhs_IntegrateByTime
         // which will scale nextFunctionEvaluation back to momentum coordinates,
         // so we have to scale back to velocity coordinates before we store it:
         for (int i = 3; i < 6; i ++)
            nextFunctionEvaluation[i] /= m_fEq -> FMass();

         // Getting number of function calls used so far:
         const G4CachedMagneticField *myField = (G4CachedMagneticField*)
                           ( BaseStepper::GetEquationOfMotion() -> GetFieldObj() );
         G4int no_function_calls = myField -> GetCountCalls();

         // Storing all this with StepTracker:
         BaseStepper::mTracker -> RecordResultOfStepper(yIn, dydx_copy,
                                                        yOutput, nextFunctionEvaluation,
                                                        hstep, // Supposed to be in arclength units.
                                                        no_function_calls);

         BaseStepper::mTracker -> UnArmTracker();
      }
   }
#endif

   // Have to convert back to momentum coordinates.
   for (int i = 3; i < 6; i ++) {
      yOutput[i] *= m_fEq -> FMass();
      yError[i] *= m_fEq -> FMass();      // Temp check to see if we have to scale last 3 coordinates of the error also.
   }
}

template <class BaseStepper>
inline
G4double MagIntegratorStepper_byTime<BaseStepper>::DistChord() const{

#ifdef TRACKING
   G4int no_function_calls_before_aux_stepper =
           (( G4CachedMagneticField* )( BaseStepper::mTracker -> getStepper() -> GetEquationOfMotion() -> GetFieldObj() ))
                                                  -> GetCountCalls();

#endif

   G4double dist_chord = BaseStepper::DistChord();

#ifdef TRACKING

   BaseStepper::mTracker -> no_function_calls_used_by_DistChord +=
          (( G4CachedMagneticField* )( BaseStepper::mTracker -> getStepper() -> GetEquationOfMotion() -> GetFieldObj() ))
             -> GetCountCalls() - no_function_calls_before_aux_stepper;
#endif

   return dist_chord;
}

template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::MagIntegratorStepper_byTime(
                                                         G4Mag_EqRhs *EquationRhs,
                                                         G4int numberOfVariables,
                                                         G4int numStateVariables)
: BaseStepper( EquationRhs, numberOfVariables )
{
   m_fEq = EquationRhs;
   for (int i = 0; i < NO_STATE_VARIABLES; i ++)
      yIn[i] = 0.;
}
template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::~MagIntegratorStepper_byTime() {
}

#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_ */

