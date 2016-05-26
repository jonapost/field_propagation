//
// class Comparator
//
// Class description:
//
// Comparator is the class for cross-checking Geant4 steppers accuracy.
// Their accuracy is checked by calling Compare() function. Parameter
// useDriver if true forces it to use G4MagInt_Driver, if false stepper's
// errors are not taken to account.
//
// History:
// - Created: D. Sorokin, May 11th 2016.
// --------------------------------------------------------------------

#ifndef COMPARATOR_HH
#define COMPARATOR_HH

#include "G4Types.hh"
#include "G4coutDestination.hh"
#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Proton.hh"
#include "G4DynamicParticle.hh"

#include "ChordFinder.hh"

enum mode{
    Default,
    Verbose,
    Silent
};

class Comparator{
public:
    Comparator(G4DynamicParticle* pDynParticle, G4MagneticField* pfield);
    ~Comparator();

    //main functions
    void CrossCheck(const G4double* const testData, const G4double* const refData, G4int mode = Default);
    template <class testStepper, class refStepper>
    void Compare(const G4double stepLen, const G4int NSteps, const bool useDriver, const G4int verbosity);

    //emulates FindNextChord()
    template <class testDriver, class refDriver>
    void CompareDrivers(const G4double stepLen, const G4int NSteps, const G4int verb);

    //setters
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


template <class testStepper, class refStepper>
void Comparator::Compare(const G4double stepLen, const G4int NSteps, const bool useDriver, const G4int verbosity){
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

//template <class testDriver, class refDriver>
//void CompareDrivers(const G4double stepLen, const G4int NSteps, const G4int verb){
//    //class T_Field, class T_Equation, class T_Stepper, class T_Driver
//    ChordFinder<G4MagneticField, G4EquationOfMotion, > chordFinder(new testDriver);




#endif
