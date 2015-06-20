/*
 * MagIntegratorStepperbyTime.hh
 *
 *  Created on: Jun 20, 2015
 *      Author: jason
 */

#ifndef MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_
#define MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_

#include "G4Mag_EqRhs.hh"
#include "G4MagIntegratorStepper.hh"

#include <iostream>
using namespace std;

template <class BaseStepper>

class MagIntegratorStepper_byTime : public G4MagIntegratorStepper {
public:
   inline
   MagIntegratorStepper_byTime(G4Mag_EqRhs *EquationRhs,
                                 G4int numberOfVariables = 6,
                                 G4int numStateVariables = 12);
   virtual ~MagIntegratorStepper_byTime();

   inline void ComputeRightHandSide(const G4double yIn[], G4double dydx[]);

   inline void Stepper(const G4double yInput[],
                       const G4double dydx[],
                             G4double hstep,
                             G4double yOutput[],
                             G4double yError [] );

   inline G4double  DistChord() const;
   inline G4int IntegratorOrder() const;


private:
   G4double yIn[10];
   G4double mass, inv_mass; // G4MagIntegratorStepper doesn't have these fields by default

   G4Mag_EqRhs * m_fEq;
   BaseStepper *baseStepper;


};

#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_ */

template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::ComputeRightHandSide(const G4double yIn[], G4double dydx[]) {
   baseStepper->ComputeRightHandSide(yIn, dydx);
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
   for (int i = 0; i < 10; i ++)
      yIn[i] = yInput[i];
   for (int i = 3; i < 6; i ++)
      yIn[i] *= inv_mass;
   baseStepper->Stepper( yInput, dydx, hstep, yOutput, yError );
   for (int i = 3; i < 6; i ++)
      yOutput[i] *= mass;
}

template <class BaseStepper>
inline
G4double MagIntegratorStepper_byTime<BaseStepper>::DistChord() const{
   return baseStepper->DistChord();
}

template <class BaseStepper>
inline
G4int MagIntegratorStepper_byTime<BaseStepper>::IntegratorOrder() const{
   return baseStepper->IntegratorOrder();
}

template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::MagIntegratorStepper_byTime(
                                                                  G4Mag_EqRhs *EquationRhs,
                                                                  G4int numberOfVariables,
                                                                  G4int numStateVariables)
: G4MagIntegratorStepper(EquationRhs, numberOfVariables, numStateVariables)
{
   m_fEq = EquationRhs;
   baseStepper = new BaseStepper(EquationRhs);  //, numberOfVariables, numStateVariables);
   for (int i = 0; i < 10; i ++)
      yIn[i] = 0.;
   mass = mass = m_fEq->FMass();
   inv_mass = 1. / mass;
}
template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::~MagIntegratorStepper_byTime() {
   delete baseStepper;
}



