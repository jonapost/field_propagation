#ifndef DetectorConstruction_HH
#define DetectorConstruction_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4RevisedChordFinder.hh"
#include "G4CachedMagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4Material.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();


    G4CachedMagneticField* GetField() const
    {return fpField;}

private:
    void ConstructField();
    void AddCube(const char* name,
                 const G4ThreeVector &halfSize,
                 const G4ThreeVector& pos,
                 G4LogicalVolume *mother,
                 G4Material *material);

private:
    G4CachedMagneticField* fpField;
    G4UniformMagField* field;
    G4RevisedChordFinder* fpChordFinder;
    G4double  fMinChordStep;   // Minimum Step for chord
};



#endif

