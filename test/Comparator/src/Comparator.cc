#include <Comparator.hh>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <G4MagneticField.hh>
using namespace CLHEP;

#define save2File

Comparator::Comparator(G4DynamicParticle* pDynParticle, G4MagneticField* pfield):
    testTrack(nullptr),
    refTrack(nullptr),
    field(pfield),
    equation(nullptr),
    dynParticle(pDynParticle),
    pos(0,0,0),
    hmin(1e-4),
    precision(1e-4),
    diffSteps(0),
    maxDiff(0)
{
    initTracks();
    initEquation();
}

Comparator::~Comparator(){
    if (testTrack != nullptr) delete testTrack;
    if (refTrack != nullptr) delete refTrack;
    if (field != nullptr) delete field;
    if (equation != nullptr) delete equation;
    if (dynParticle != nullptr) delete dynParticle;
}

void Comparator::CrossCheck(const G4double * const testData, const G4double * const refData, G4int mode){
#ifdef save2File
    static std::ofstream out("out.txt");
    out << testData[0] <<"  "<< testData[1] <<"  "<< testData[2] << G4endl;
#endif
    G4double h2 = 0, epsh = 0, p2 = 0, epsp = 0;
    for (G4int i = 0; i < 3; ++i){
        epsh += sqr(testData[i] - refData[i]);
        h2 += sqr(refData[i]);
        p2 += sqr(refData[i+3]);
        epsp += sqr(testData[i+3] - refData[i+3]);
    }
    epsh /= h2;
    epsp /= p2;
    G4double eps = std::max(epsh,epsp)/(precision*precision);
    if (mode == Default){
        if (eps > 1)
            G4cout << eps << G4endl;
    }
    else if (mode == Verbose){
        if (eps > 1){
            G4cout<<"("<<testData[0]<<","<<testData[1]<<","<<testData[2]<<")   "
                  <<"("<<refData[0] <<","<<refData[1] <<","<<refData[2]<<")   "
                  << G4endl;
        }
    }
    else if (mode == Silent){
        //accamulate statistics
        if (eps > 1){
            ++diffSteps;
            maxDiff = std::max(maxDiff,eps);
        }
    }
    else{
        G4cout<<"wrong mode <"<<mode<<"> \n";
    }
}

void Comparator::initTracks(){

    if (testTrack != nullptr) delete testTrack;
    if (refTrack != nullptr) delete refTrack;

    testTrack = new  G4FieldTrack(pos,
                              0,
                              dynParticle->GetMomentumDirection(),
                              dynParticle->GetKineticEnergy(),
                              dynParticle->GetMass(),
                              dynParticle->GetCharge(),
                              dynParticle->GetPolarization());
    refTrack = new G4FieldTrack(*testTrack);
}

void Comparator::initEquation(){
    if (equation != nullptr) delete equation;

    equation = new G4Mag_UsualEqRhs(field);
    const G4ChargeState chargeState(dynParticle->GetCharge(), dynParticle->GetSpin(), dynParticle->GetMagneticMoment());
    equation->SetChargeMomentumMass(chargeState, dynParticle->GetMomentum().mag(), dynParticle->GetMass());
}

void Comparator::setParticle(G4DynamicParticle * pDynParticle){
    if (pDynParticle != dynParticle){
        delete dynParticle;
        dynParticle = pDynParticle;
        initTracks();
        initEquation();
    }
}

void Comparator::setField(G4MagneticField * pfield){
    if (pfield != field){
        delete field;
        field = pfield;
        initEquation();
    }
}

void Comparator::setPostition(const G4ThreeVector &newPos){
    pos = newPos;
    initTracks();
}

void Comparator::setMinDriverStep(const G4double newHmin){
    hmin = newHmin;
}
