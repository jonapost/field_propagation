#include "BSStepper.hh"

BSStepper::BSStepper(G4EquationOfMotion *Equation,
                     G4int  numIntegrationVariables,
                     G4int  numStateVariables):
    G4MagIntegratorStepper(Equation,numIntegrationVariables,numStateVariables)
{

}

BSStepper::~BSStepper(){

}


void  BSStepper::Stepper(const G4double /*y*/[],
                         const G4double /*dydx*/[],
                         G4double /*h*/,
                         G4double /*yout*/[],
                         G4double /*yerr*/[]  ){
    G4cout<<"BSStepper::Stepper should never been called! \n";
    throw("wrong call!");
}


G4double  BSStepper::DistChord() const {
    G4cout<<"BSStepper::DistChord should never been called! \n";
    throw("wrong call!");
    return -1;
}



G4int BSStepper::IntegratorOrder() const {
    G4cout<<"BSStepper::IntegratorOrder should never been called! \n";
    throw("wrong call!");
    return -1;
}
