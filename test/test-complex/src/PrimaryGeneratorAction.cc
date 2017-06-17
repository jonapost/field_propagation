//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file persistency/gdml/G01/src/G01PrimaryGeneratorAction.cc
/// \brief Implementation of the G01PrimaryGeneratorAction class
//
//
// $Id: G01PrimaryGeneratorAction.cc 68025 2013-03-13 13:43:46Z gcosmo $
//
//
// Original Geant4 implementation was:
//   examples/extended/persistency/gdml/G01/src/G01PrimaryGeneratorAction.cc 

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"


#ifdef USE_VECGEOM_NAVIGATOR
#undef USE_VECGEOM_NAVIGATOR
#define RESTORE_USE_VECGEOM_NAVIGATOR
#endif

#include "HepMCGenerator.h"
#include "PrimaryGeneratorAction.hh"

#ifdef RESTORE_USE_VECGEOM_NAVIGATOR
#define USE_VECGEOM_NAVIGATOR
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction() 
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   fHepMCGenerator(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fParticleGun->SetParticleDefinition(
               particleTable->FindParticle(particleName="geantino"));
  fParticleGun->SetParticleEnergy(1.0*GeV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.366847304086928, 0.262356272922734, 0.892520163101229));
}

PrimaryGeneratorAction::PrimaryGeneratorAction(std::string& filename)
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   fHepMCGenerator(0)
{
 fHepMCGenerator = new HepMCGenerator(filename); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if(fParticleGun)
  delete fParticleGun;
  if(fHepMCGenerator)
  delete fHepMCGenerator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent);
  G4cout<<"-------- Primary has been generated "<<G4endl;
}
*/

// Take one HepMCGenerator event and fill Geant4 event with the primaries 
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   G4ParticleTable  *particleTable = G4ParticleTable::GetParticleTable();
   G4double         tpx, tpy, tpz, te, x0, y0, z0;
   G4int            pdg;
   
   GeantEventInfo eventinfo = fHepMCGenerator->NextEvent();
   G4PrimaryVertex  *pvertex = new G4PrimaryVertex(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert, eventinfo.tvert);

   G4int ntracks = eventinfo.ntracks;
   for(G4int ip = 0; ip < ntracks; ++ip){
     fHepMCGenerator->GetTrack(ip, tpx, tpy, tpz, te, x0, y0, z0, pdg);
     G4ParticleDefinition *particleDef = particleTable->FindParticle(pdg);
     pvertex->SetPrimary(new G4PrimaryParticle(particleDef, tpx, tpy, tpz, te));
   }
   anEvent->AddPrimaryVertex(pvertex);
}

