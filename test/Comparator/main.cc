#include "Comparator.hh"

#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4CachedMagneticField.hh"

//first order explicit RK methods
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4ExplicitEuler.hh"

//second order explicit RK methods
#include "G4BogackiShampine23.hh"

//3rd order explicit RK methods
#include "G4RKG3_Stepper.hh"
#include "G4DoLoMcPriRK34.hh"

//4th order explicit RK methods
#include "G4ClassicalRK4.hh"
#include "G4ConstRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4TsitourasRK45.hh"
#include "G4DormandPrince745.hh"
#include "G4BogackiShampine45.hh"

//5th order explicit RK methods
#include "G4DormandPrinceRK56.hh"

//7th order explicit RK methods
#include "G4DormandPrinceRK78.hh"

// exact steppers
#include "G4ExactHelixStepper.hh"


#include <iostream>
#include <memory>

using namespace CLHEP;


int main()
{
    auto dynParticle =
        std::make_unique<G4DynamicParticle>(
            G4Proton::Definition(), G4ThreeVector(1, 0.01, 0.01), 1*GeV);

    G4cout << "cat " << G4endl;
    auto field = std::make_unique<G4UniformMagField>(G4ThreeVector(0, 0, 1*tesla));
    //auto field = std::make_unique<G4QuadrupoleMagField>(0.1 * tesla / meter);
    auto magneticField = std::make_shared<G4CachedMagneticField>(field.get(), 0);

    Comparator comparator(std::move(dynParticle), magneticField);

    //comparator.compare<G4ExactHelixStepper, G4ExactHelixStepper>(
    //   100*cm, 10000, Comparator::Mode::SaveError);

    comparator.setPrecision(1);

    comparator.compareDriver<StandardDriver, G4MagInt_Driver, G4ClassicalRK4>(
       10000*cm,Comparator::Mode::SaveTrack);

    //comparator.CompareWithBS<G4CashKarpRKF45>(1000*m, Comparator::Mode::Verbose);

    return 0;
}
