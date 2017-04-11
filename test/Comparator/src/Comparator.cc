#include "Comparator.hh"

#include "G4MagneticField.hh"

#include <cmath>
#include <iomanip>
#include <fstream>

using namespace CLHEP;

Comparator::Comparator (std::unique_ptr<G4DynamicParticle> dynParticle,
    std::shared_ptr<G4CachedMagneticField> field):
    ftestTrack(nullptr),
    frefTrack(nullptr),
    fTestField(field),
    fRefField(std::make_shared<G4CachedMagneticField>(*fTestField)),
    fTestEquation(nullptr),
    fRefEquation(nullptr),
    fdynParticle(std::move(dynParticle)),
    fstartPosition(0,0,0),
    fhmin(1e-4),
    fprecision(1e-4),
    fdiffSteps(0),
    fmaxDiff(0),
    ffout("out.txt")
{
    initialize();
}

Comparator::~Comparator()
{
    ffout.close();
}

void Comparator::initialize()
{
    ftestTrack = std::make_unique<G4FieldTrack>(fstartPosition, 0,
        fdynParticle->GetMomentumDirection(), fdynParticle->GetKineticEnergy(),
            fdynParticle->GetMass(), fdynParticle->GetCharge(),
                fdynParticle->GetPolarization());

    frefTrack = std::make_unique<G4FieldTrack>(*ftestTrack);

    fTestEquation = std::make_shared<G4Mag_UsualEqRhs>(fTestField.get());
    fRefEquation = std::make_shared<G4Mag_UsualEqRhs>(fRefField.get());

    const G4ChargeState chargeState(fdynParticle->GetCharge(),
        fdynParticle->GetSpin(), fdynParticle->GetMagneticMoment());

    fTestEquation->SetChargeMomentumMass(chargeState,
        fdynParticle->GetMomentum().mag(), fdynParticle->GetMass());

    fRefEquation->SetChargeMomentumMass(chargeState,
        fdynParticle->GetMomentum().mag(), fdynParticle->GetMass());
}

/*
 * data[0..2] - position
 * data[3..5] - momentum
 *
 * */

void Comparator::crossCheck(const G4double* const testData,
    const G4double* const refData, Mode mode)
{
    G4double h2 = 0, epsh = 0, p2 = 0, epsp = 0;
    for (G4int i = 0; i < 3; ++i) {
        epsh += sqr(testData[i] - refData[i]);
        h2 += sqr(refData[i]);
        p2 += sqr(refData[i+3]);
        epsp += sqr(testData[i+3] - refData[i+3]);
    }
    epsh /= h2;
    epsp /= p2;
    G4double eps = std::max(epsh,epsp)/sqr(fprecision);

    if (eps > 1) {
        if (mode == Mode::Default) {
            G4cout << eps << G4endl;
        } else if (mode == Mode::Verbose) {
            G4cout << G4ThreeVector (testData[0], testData[1], testData[2]) << "  "
                << G4ThreeVector (refData[0], refData[1], refData[2]) << G4endl;
        } else if (mode == Mode::Silent) {
            ++fdiffSteps;
            fmaxDiff = std::max(fmaxDiff, eps);
        }
    }

    if (mode == Mode::SaveTrack) {
        ffout << testData[0] <<"  "<< testData[1] <<"  "<< testData[2] << G4endl;
    }
    if (mode == Mode::SaveError) {
        ffout << eps << G4endl;
    }
}

void Comparator::setParticle(std::unique_ptr<G4DynamicParticle> dynParticle)
{
    if (fdynParticle != dynParticle) {
        fdynParticle = std::move(dynParticle);
        initialize();
    }
}

void Comparator::setField(std::shared_ptr<G4CachedMagneticField> field)
{
    if (fTestField != field) {
        fTestField = field;
        fRefField = std::make_shared<G4CachedMagneticField>(*fTestField);
        initialize();
    }
}

void Comparator::setStartPostition(const G4ThreeVector& position)
{
    fstartPosition = position;
    initialize();
}

void Comparator::setMinDriverStep(const G4double hmin)
{
    fhmin = hmin;
}

void Comparator::setPrecision(const G4double precision)
{
    fprecision = precision;
}
