#include "Comparator.hh"

#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
#include "CMSMagneticField.hh"
#include "G4CachedMagneticField.hh"

//first order explicit RK methods
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4ExplicitEuler.hh"

//second order explicit RK methods
#include "G4BogackiShampine23.hh"

//3rd order explicit RK methods
#include "G4RKG3_Stepper.hh"

//4th order explicit RK methods
//#include "G4ClassicalRK4.hh"
#include "G4ConstRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4DormandPrince745.hh"
#include "G4BogackiShampine45.hh"

//5th order explicit RK methods

//7th order explicit RK methods

// exact steppers
#include "G4ExactHelixStepper.hh"


#include "RK547FEq1.hh"

#include <iostream>
#include <memory>

using namespace CLHEP;


int main()
{
    auto dynParticle =
        std::make_unique<G4DynamicParticle>(
            G4Proton::Definition(), G4ThreeVector(1, 0, 2).unit(), 0.01*GeV);

    //auto field = std::make_unique<G4UniformMagField>(G4ThreeVector(0, 0, 1*tesla));
    //auto field = std::make_unique<G4QuadrupoleMagField>(1 * tesla / meter);
    auto field = std::make_unique<CMSMagneticField>("cmsmagfield2015.txt");
    auto magneticField = std::make_shared<G4CachedMagneticField>(field.get(), 0);

    Comparator comparator(std::move(dynParticle), magneticField);

   // comparator.compare<G4ClassicalRK4, /*G4ExactHelixStepper*/G4ClassicalRK4>(
   //    100*cm, 10000, Comparator::Mode::SaveTrack);

    comparator.setPrecision(1);
    comparator.setStartPostition(G4ThreeVector(100, 100, -16000));

    comparator.compareDriver<GustafssonDriver, G4MagInt_Driver, G4DormandPrince745>(
       40000*mm, Comparator::Mode::SaveTrack);

    //comparator.CompareWithBS<G4CashKarpRKF45>(1000*m, Comparator::Mode::Verbose);
    return 0;
}
