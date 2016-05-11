#ifndef COMPARATOR_HH
#define COMPARATOR_HH

#include <array>
#include <cstring>
#include <fstream>
#include <memory>
#include <iomanip>

#include "G4Types.hh"
#include "G4coutDestination.hh"
#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Proton.hh"
#include "G4DynamicParticle.hh"

enum mode{
    Default,
    Verbose,
    Silent
};

class Comparator{
public:
    Comparator(G4DynamicParticle* pDynParticle, G4MagneticField* pfield):
        testTrack(nullptr),
        refTrack(nullptr),
        field(pfield),
        equation(nullptr),
        dynParticle(pDynParticle),
        pos(0,0,0),
        hmin(1e-6),
        precision(1e-4),
        diffSteps(0),
        maxDiff(0)
    {
        initTracks();
        initEquation();
    }
    ~Comparator();

    void CrossCheck(const G4double* const testData, const G4double* const refData, G4int mode = Default);
    template <class testStepper, class refStepper>
    void Compare(const G4double stepLen, const G4int NSteps, const bool useDriver, const G4int verbosity){
        G4double yIn[N],yOut[N],yErr[N],dydx[N];
        G4double yInRef[N],yOutRef[N],yErrRef[N],dydxRef[N];

        testStepper testSt(equation);
        refStepper refSt(equation);
        diffSteps = 0;
        if (verbosity == Verbose){
            G4cout<<"test position                ref position \n ";
        }

        if (useDriver){
            G4MagInt_Driver testDr(hmin,&testSt);
            G4MagInt_Driver refDr(hmin,&refSt);
            for (G4int i = 0; i < NSteps; ++i){
                testDr.AccurateAdvance(*testTrack,stepLen,precision,stepLen);
                refDr.AccurateAdvance(*refTrack,stepLen,precision,stepLen);

                testTrack->DumpToArray(yOut);
                refTrack->DumpToArray(yOutRef);

                CrossCheck(yOut,yOutRef,verbosity);
            }
        }
        else{
            for (G4int i = 0; i < NSteps; ++i){
                testTrack->DumpToArray(yIn);
                refTrack->DumpToArray(yInRef);

                testSt.RightHandSide(yIn, dydx);
                testSt.Stepper(yIn,dydx,stepLen,yOut,yErr);

                refSt.RightHandSide(yInRef, dydxRef);
                refSt.Stepper(yInRef,dydxRef,stepLen,yOutRef,yErrRef);

                testTrack->LoadFromArray(yOut,N);
                refTrack->LoadFromArray(yOutRef,N);

                CrossCheck(yOut,yOutRef,verbosity);
            }
        }

        if (verbosity == Silent){
            G4cout<<"diffSteps# "<<diffSteps<<" maxDiff "<<maxDiff<<G4endl;
        }
    }

    void setParticle(G4DynamicParticle* pDynParticle);
    void setField(G4MagneticField* pfield);
    void setPostition(const G4ThreeVector& newPos);
    void setMinDriverStep(const G4double newHmin);
    void setPrecision(const G4double newPrec);



private:

    void initTracks();
    void initEquation();


    G4FieldTrack* testTrack;
    G4FieldTrack* refTrack;
    G4MagneticField *field;
    G4Mag_UsualEqRhs *equation;

    //
    G4DynamicParticle* dynParticle;
    G4ThreeVector pos;
    G4double hmin;
    G4double precision;

    static const G4int N = G4FieldTrack::ncompSVEC;

    G4int diffSteps;
    G4double maxDiff;
};

#endif
