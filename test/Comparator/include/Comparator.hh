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

#pragma once

#include "G4Types.hh"
#include "G4coutDestination.hh"
#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4CachedMagneticField.hh"

#include <fstream>
#include <memory>

class Comparator {
public:
    enum class  Mode {
        Default,
        Verbose,
        Silent,

        SaveTrack,
        SaveError
    };

    Comparator(std::unique_ptr<G4DynamicParticle> dynParticle,
        std::shared_ptr<G4CachedMagneticField> field);

    ~Comparator();

    //main functions
    void crossCheck(const G4double* const testData,
        const G4double* const refData, Mode mode = Mode::Default);

    template <class TestStepper, class RefStepper>
    void compare(
        const G4double stepLen, const G4int nSteps, Mode mode);
    
    template <class TestDriver, class RefDriver, class Stepper>
    void compareDriver(G4double length, Mode mode);

    //compare with Bulirsch-Stoer driver
    template <class RefStepper>
    void CompareWithBS(const G4double path, const G4int verb);

    //setters
    void setParticle (std::unique_ptr<G4DynamicParticle> dynParticle);
    void setField (std::shared_ptr<G4CachedMagneticField> field);
    void setStartPostition (const G4ThreeVector& position);
    void setMinDriverStep (const G4double hmin);
    void setPrecision (const G4double precision);

private:
    void initialize ();

    std::unique_ptr<G4FieldTrack> ftestTrack;
    std::unique_ptr<G4FieldTrack> frefTrack;

    std::shared_ptr<G4CachedMagneticField> fTestField;
    std::shared_ptr<G4CachedMagneticField> fRefField;
    std::shared_ptr<G4Mag_UsualEqRhs> fTestEquation;
    std::shared_ptr<G4Mag_UsualEqRhs> fRefEquation;

    std::unique_ptr<G4DynamicParticle> fdynParticle;
    G4ThreeVector fstartPosition;
    G4double fhmin;
    G4double fprecision;

    G4int fdiffSteps;
    G4double fmaxDiff;

    std::ofstream ffout;

    static const G4int fncomp = G4FieldTrack::ncompSVEC;
};


#include "Comparator.icc"
