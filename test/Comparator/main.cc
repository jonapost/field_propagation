#include <iostream>

#include "Comparator.hh"
#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"

#include "G4ClassicalRK4.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ExactHelixStepper.hh"
#include "G4CashKarpRKF45.hh"
#include "G4ExplicitEuler.hh"


using namespace std;
using namespace CLHEP;




int main(){
    G4DynamicParticle* dynParticle = new G4DynamicParticle(G4Proton::Definition(),G4ThreeVector(1,0.01,0.01),1*GeV);
    G4MagneticField* field = new G4UniformMagField(G4ThreeVector(0, 0, 1*tesla));
    //G4MagneticField* field = new G4QuadrupoleMagField(0.01*tesla/meter);
    Comparator comparator(dynParticle,field);
    //comparator.Compare<G4ExplicitEuler,G4CashKarpRKF45>(10*cm,10000,false,Verbose);
    comparator.CompareWithBS<G4CashKarpRKF45>(100*m,Verbose);

    return 0;
}
