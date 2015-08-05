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
   G4double yIn[NO_STATE_VARIABLES],
            dydx_copy[NO_STATE_VARIABLES],
            cached_dydx[NO_STATE_VARIABLES],
            last_function_evaluation[NO_STATE_VARIABLES];

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

#ifdef TRACKING
   if ( ! BaseStepper::mTracker -> get_within_AdvanceChordLimited() )
      cout << "In Stepper, but not in AdvanceChordLimited()" << endl;
#endif

#ifdef BACK_TRACKING

   std::vector<G4double> *last_pos_vector, *second_to_last_pos_vector;
   G4int buffer_length;

   while ( ! BaseStepper::mTracker ->
         check_that_wasnt_disgarded_by_Propagator( yInput ) ) {

      buffer_length = BaseStepper::mTracker -> getBufferLength();
      // Because of current bug with BogackiShampine45 (and possibly others):
      assert( buffer_length >= 2);

      BaseStepper::mTracker -> set_last_time_val_was_accepted( true );

      last_pos_vector =
            &( BaseStepper::mTracker -> get_buffer_ptr() -> at( buffer_length - 1 ) );
      second_to_last_pos_vector =
            &( BaseStepper::mTracker -> get_buffer_ptr() -> at( buffer_length - 2 ) );

      if (     ( G4ThreeVector(yInput[0], yInput[1], yInput[2])
                 - G4ThreeVector(second_to_last_pos_vector -> at(POSITION_SLOT + 0),
                                 second_to_last_pos_vector -> at(POSITION_SLOT + 1),
                                 second_to_last_pos_vector -> at(POSITION_SLOT + 2)) ).mag()
            >=
               ( G4ThreeVector(yInput[0], yInput[1], yInput[2])
                 - G4ThreeVector(last_pos_vector -> at(POSITION_SLOT + 0),
                                 last_pos_vector -> at(POSITION_SLOT + 1),
                                 last_pos_vector -> at(POSITION_SLOT + 2)) ).mag()     )
      {
         break; // We are done searching backwards for a good match.
      }

      BaseStepper::mTracker -> get_buffer_ptr() -> pop_back();
      BaseStepper::mTracker -> get_no_function_calls_buffer() -> pop_back();

   }
#endif

   // Have to copy because yInput and dydx are constant in the function signature...
   for (int i = 0; i < 10; i ++){
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
         BaseStepper::mTracker -> RecordResultOfStepper(yOutput,
                                     nextFunctionEvaluation, no_function_calls);

         BaseStepper::mTracker -> UnArmTracker();
      }
   }
#endif

   // Have to convert back to momentum coordinates.
   for (int i = 3; i < 6; i ++)
      yOutput[i] *= m_fEq -> FMass();

}

template <class BaseStepper>
inline
G4double MagIntegratorStepper_byTime<BaseStepper>::DistChord() const{

   return BaseStepper::DistChord();
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

