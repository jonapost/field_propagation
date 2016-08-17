#ifndef SteppingAction_HH
#define SteppingAction_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "EventAction.hh"

#include "G4CachedMagneticField.hh"
#include "G4ParticleGun.hh"

class G4LogicalVolume;

/// Stepping action class
/// 

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();
    virtual ~SteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*step);

  private:

    G4CachedMagneticField* fpField;
    G4int ncalls;

    G4ParticleGun* fParticleGun;

    G4double time;
};



#endif
