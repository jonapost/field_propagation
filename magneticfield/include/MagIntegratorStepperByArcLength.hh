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
#ifndef MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYARCLENGTH_HH_
#define MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYARCLENGTH_HH_


#include "G4Mag_EqRhs.hh"
#include "G4MagIntegratorStepper.hh"

#include <iostream>
#include <assert.h>
using namespace std;

#define NO_STATE_VARIABLES 12


// G4MagIntegratorStepper_byArcLength, (a wrapper class so that CashKarp, ClassicalRK4 and others
// can interface with StepTracker when Stepper() is called, and to store total time).

template <class BaseStepper>

class MagIntegratorStepper_byArcLength : public BaseStepper {
public:

   MagIntegratorStepper_byArcLength(G4Mag_EqRhs *EquationRhs,
                                    G4int numberOfVariables = 6,
                                    G4int numStateVariables = 12);
   virtual ~MagIntegratorStepper_byArcLength();

   inline void ComputeRightHandSide(const G4double yInput[], G4double dydx[]);

   inline void Stepper(const G4double yInput[],
                       const G4double dydx[],
                             G4double hstep,
                             G4double yOutput[],
                             G4double yError [] );

   inline G4double DistChord() const;

private:
   // Used simply for parsing for storage to StepTracker (because Stepper()
   // takes yInput[] and dydx[] as const input.
   G4double yIn[NO_STATE_VARIABLES],
            dydx_copy[NO_STATE_VARIABLES];

   G4double nextFunctionEvaluation[NO_STATE_VARIABLES]; // 8 for safety (only need 6).

   G4Mag_EqRhs  *m_fEq;
};


template <class BaseStepper>
inline
void MagIntegratorStepper_byArcLength<BaseStepper>::ComputeRightHandSide(const G4double yInput[], G4double dydx[]) {

   BaseStepper::ComputeRightHandSide(yInput, dydx);

}

template <class BaseStepper>
inline
void MagIntegratorStepper_byArcLength<BaseStepper>::Stepper(const G4double yInput[],
      const G4double dydx[],
            G4double hstep,
            G4double yOutput[],
            G4double yError [] ) {

   BaseStepper::Stepper( yInput, dydx, hstep, yOutput, yError );
   //( dynamic_cast<BaseStepper*>( this )) -> BaseStepper::Stepper( yIn, dydx_copy, hstep, yOutput, yError );

#ifdef TRACKING
   if ( BaseStepper::mTracker -> get_within_AdvanceChordLimited() ) {
      // Within AdvancedChordLimited, so we want to record.
      if ( BaseStepper::mTracker -> isArmed() ) {
         // BaseStepper::mTracker -> RecordResultOfStepper(yIn, dydx_copy);

         // nextFunctionEvaluation is done at the right endpoint of the integration step.
         // If this works should make nextFunctionEvaluation a member variable:

         BaseStepper::ComputeRightHandSide(yOutput, nextFunctionEvaluation);

         // Want to store velocity, not momentum.
         for (int i = 3; i < 6; i ++)
                     yOutput[i] /= m_fEq -> FMass();
         // Divide by mass to get acceleration:
         for (int i = 3; i < 6; i ++)
            nextFunctionEvaluation[i] /= m_fEq -> FMass();

         // Since we are also storing the values at the beginning of the step,
         // and because we have to convert to velocity coordinates to store in StepTracker:
         for (int i = 0; i < 6; i ++){
            yIn[i] = yInput[i];
            dydx_copy[i] = dydx[i];
         }
         for (int i = 3; i < 6; i ++){
            yIn[i] *= 1. / m_fEq -> FMass();
            dydx_copy[i] *= 1. / m_fEq -> FMass();
         }

         const G4CachedMagneticField *myField = (G4CachedMagneticField*)
                                       ( BaseStepper::GetEquationOfMotion() -> GetFieldObj() );
         G4int no_function_calls = myField -> GetCountCalls();

         BaseStepper::mTracker -> RecordResultOfStepper(yIn, dydx_copy,
                                                        yOutput, nextFunctionEvaluation,
                                                        hstep, no_function_calls); // Store as velocity (instead of mom.)

         // Change back to momentum coordinates for next round of stepper:
         for (int i = 3; i < 6; i ++)
               yOutput[i] *= m_fEq -> FMass();

         BaseStepper::mTracker -> UnArmTracker();
      }
   }
#endif

}

template <class BaseStepper>
inline
G4double MagIntegratorStepper_byArcLength<BaseStepper>::DistChord() const{
   return BaseStepper::DistChord(); // If BaseStepper uses Richardson
                                    // extrapolation (auxStepper) then this might
                                    // not be producing the right results.
}

template <class BaseStepper>
inline MagIntegratorStepper_byArcLength<BaseStepper>::MagIntegratorStepper_byArcLength(
                                                                  G4Mag_EqRhs *EquationRhs,
                                                                  G4int numberOfVariables,
                                                                  G4int numStateVariables)
: BaseStepper( EquationRhs, numberOfVariables )//, numStateVariables)
{

   // Not really needed? (It is needed for MagIntegratorStepper_byTime because we need FMass(),
   // which is not in G4EquationOfMotion).
   m_fEq = EquationRhs;
}

template <class BaseStepper>
MagIntegratorStepper_byArcLength<BaseStepper>::~MagIntegratorStepper_byArcLength() { }


#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYARCLENGTH_HH_ */
