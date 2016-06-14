// class BSStepper
//
// Class description:
//
// dummy class for BulirschStoer driver
//
// History:
// - created Dmitry Sorokin
// --------------------------------------------------------------------

#ifndef BSSTEPPER_HH
#define BSSTEPPER_HH
#include "G4MagIntegratorStepper.hh"

class BSStepper: public G4MagIntegratorStepper{
public:
    BSStepper(G4EquationOfMotion *Equation,
              G4int  numIntegrationVariables = 6,
              G4int  numStateVariables = 12);

   virtual ~BSStepper();


   virtual  void  Stepper(const G4double y[],
                          const G4double dydx[],
                          G4double h,
                          G4double yout[],
                          G4double yerr[]  );


   virtual  G4double  DistChord() const;


   virtual G4int IntegratorOrder() const;

};

#endif
