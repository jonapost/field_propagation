#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4Proton.hh"

using namespace CLHEP;

PrimaryGeneratorAction::PrimaryGeneratorAction():
    G4VUserPrimaryGeneratorAction(),fParticleGun(nullptr)
{
  G4int n_particle = 1;
  G4double energy = 1* MeV;
  G4double BField = 1*tesla;

  G4ParticleDefinition* particle = G4Proton::Definition();

  fParticleGun  = new G4ParticleGun(n_particle);

  fParticleGun->SetParticleDefinition(particle);

  G4double mass = particle->GetPDGMass();
  G4double charge = particle->GetPDGCharge();
  G4double mom2 = energy*(energy + 2*mass);
  G4ThreeVector momDir = G4ThreeVector(1,1,0).unit();

  G4double momXZ = sqrt(mom2)*sqrt(sqr(momDir.x()) + sqr(momDir.z()))*c_light;
  G4double radius = momXZ/(charge*BField*c_squared);

  fParticleGun->SetParticleMomentumDirection(momDir);
  fParticleGun->SetParticleEnergy(energy);
  fParticleGun->SetParticlePosition(G4ThreeVector(0*cm,-4.1*m,-radius));
}



PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}



void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent);
}



