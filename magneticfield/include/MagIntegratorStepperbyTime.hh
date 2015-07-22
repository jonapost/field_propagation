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

#include <iostream>
using namespace std;

//G4MagIntegratorStepper

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

   //inline void SetTrue_last_step_succeeded();

   //inline void Reset_last_step_succeeded();

   inline G4double DistChord() const;
   //inline G4int IntegratorOrder() const;


   //G4double mass, inv_mass;   // G4MagIntegratorStepper doesn't have these fields by default
                              // only public for quickness!
   //BaseStepper *baseStepper;
   //G4MagIntegratorStepper *baseStepper;

private:
   G4double yIn[10], dydx_copy[10], cached_dydx[10], last_function_evaluation[10];
   //G4double mass, inv_mass; // G4MagIntegratorStepper doesn't have these fields by default

   // G4bool last_step_succeeded;
   // last_step_succeeded now inherits from G4MagIntegratorStepper

   G4Mag_EqRhs  *m_fEq;

};


template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::ComputeRightHandSide(const G4double yInput[], G4double dydx[]) {
   /*if ( BaseStepper::get_last_step_succeeded() ) {
      dydx = cached_dydx;
      return;
   }
   */

   for (int i = 0; i < 10; i ++)
      yIn[i] = yInput[i];
   for (int i = 3; i < 6; i ++)
      yIn[i] *= 1. / m_fEq -> FMass();

   BaseStepper::ComputeRightHandSide(yIn, dydx);

   // Cache (for possible FSAL):
   //for (int i = 0; i < 6; i ++)
   //   last_function_evaluation[i] = dydx[i];

   for (int i = 3; i < 6; i ++)
      dydx[i] *= m_fEq -> FMass();
   // Always feed this template class a Mag_UsualEqRhs_IntegrateByTime as EquationRhs
   //for (int i = 3; i < 6; i ++)
   //   dydx[i] *= inv_mass;
}

template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::Stepper(const G4double yInput[],
      const G4double dydx[],
            G4double hstep,
            G4double yOutput[],
            G4double yError [] ) {

   // Have to copy because yInput and dydx are constant in the function signature.
   for (int i = 0; i < 10; i ++){
      yIn[i] = yInput[i];
      dydx_copy[i] = dydx[i];
   }
   for (int i = 3; i < 6; i ++){
      yIn[i] *= 1. / m_fEq -> FMass();
      dydx_copy[i] *= 1. / m_fEq -> FMass();
   }

   BaseStepper::Stepper( yIn, dydx_copy, hstep, yOutput, yError );

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

   for (int i = 3; i < 6; i ++)
      yOutput[i] *= m_fEq -> FMass();

   //yOutput[7] = yInput[7] + hstep;

}

template <class BaseStepper>
inline
G4double MagIntegratorStepper_byTime<BaseStepper>::DistChord() const{
   return BaseStepper::DistChord(); // If BaseStepper uses Richardson
                                    // extrapolation (auxStepper) then this might
                                    // not be producing the right results.
}


template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::MagIntegratorStepper_byTime(
                                                                  G4Mag_EqRhs *EquationRhs,
                                                                  G4int numberOfVariables,
                                                                  G4int numStateVariables)
: BaseStepper( EquationRhs, numberOfVariables )//, numStateVariables)
{

   //m_fEq = EquationRhs;
   // baseStepper = new BaseStepper(EquationRhs);  //, numberOfVariables, numStateVariables);
   for (int i = 0; i < 10; i ++)
      yIn[i] = 0.;

}
template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::~MagIntegratorStepper_byTime() {
   // delete baseStepper;

}


#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_ */

