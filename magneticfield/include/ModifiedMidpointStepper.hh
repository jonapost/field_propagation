// class ModifiedMidpointStepper
//
// Class description:
//
// dummy class for BulirschStoer driver
//
// History:
// - created Dmitry Sorokin
// --------------------------------------------------------------------

#ifndef ModifiedMidpointStepper_HH
#define ModifiedMidpointStepper_HH
#include "G4MagIntegratorStepper.hh"
#include "ModifiedMidpointDenseOut.hh"
#include "G4FieldTrack.hh"

class ModifiedMidpointStepper: public G4MagIntegratorStepper{
public:
    ModifiedMidpointStepper(G4EquationOfMotion *Equation,
              G4int  numIntegrationVariables = 6,
              G4int  numStateVariables = 12);

   virtual ~ModifiedMidpointStepper();


   virtual  void  Stepper(const G4double y[],
                          const G4double dydx[],
                          G4double hstep,
                          G4double yout[],
                          G4double yError[]  );


   virtual  G4double  DistChord() const;


   virtual G4int IntegratorOrder() const;

private:
    ModifiedMidpointDenseOut denseMidpoint;

    const int interval_sequence[2];
    const G4double fcoeff;

    G4double yIn[G4FieldTrack::ncompSVEC],
             yOut[2][G4FieldTrack::ncompSVEC],
             yMid[2][G4FieldTrack::ncompSVEC],
             derivs[2][6][G4FieldTrack::ncompSVEC];

};

#endif
