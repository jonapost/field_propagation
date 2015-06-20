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

   BaseStepper *baseStepper;


};

#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_ */

template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::ComputeRightHandSide(const G4double yIn[], G4double dydx[]) {
   baseStepper->ComputeRightHandSide(yIn, dydx);
}

template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::Stepper(const G4double yInput[],
      const G4double dydx[],
            G4double hstep,
            G4double yOutput[],
            G4double yError [] ) {
   baseStepper->Stepper( yInput, dydx, hstep, yOutput, yError );
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
   baseStepper = new BaseStepper(EquationRhs);  //, numberOfVariables, numStateVariables);
   // TODO Auto-generated constructor stub
}
template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::~MagIntegratorStepper_byTime() {
   delete baseStepper;
   // TODO Auto-generated constructor stub

}



