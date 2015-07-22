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
using namespace std;

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
   //G4double yIn[10], dydx_copy[10], cached_dydx[10], last_function_evaluation[10];
   //G4double mass, inv_mass; // G4MagIntegratorStepper doesn't have these fields by default

   G4Mag_EqRhs  *m_fEq;

};


template <class BaseStepper>
inline
void MagIntegratorStepper_byArcLength<BaseStepper>::ComputeRightHandSide(const G4double yInput[], G4double dydx[]) {

   BaseStepper::ComputeRightHandSide(yInput, dydx);

   // Cache:
   // for (int i = 0; i < 6; i ++)
   //   last_function_evaluation[i] = dydx[i];

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
   // BaseStepper::mTracker -> RecordResultOfStepper(yIn, dydx_copy);

   // nextFunctionEvaluation is done at the right endpoint of the integration step.
   // If this works should make nextFunctionEvaluation a member variable:
   G4double nextFunctionEvaluation[8]; // 8 for safety (only need 6).

   // Want to store velocity, not momentum, so wait to multiply yOutput[3..5] by FMass()
   // until after compute next function evaluation.
   BaseStepper::ComputeRightHandSide(yOutput, nextFunctionEvaluation);
   BaseStepper::mTracker -> RecordResultOfStepper(yOutput, nextFunctionEvaluation);
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
