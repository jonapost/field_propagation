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
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "SimplePhysicsList.hh"

#include "G4ProcessManager.hh"

#include "G4UserSpecialCuts.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4SystemOfUnits.hh"


#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4PhysicsListHelper.hh"

#ifdef USE_VECGEOM_NAVIGATOR
#undef USE_VECGEOM_NAVIGATOR
#define RESTORE_USE_VECGEOM_NAVIGATOR
#endif

#include "TotalPhysicsProcess.hh"
#include "TabulatedDataManager.hh"
#include "TPartIndex.h"

#ifdef RESTORE_USE_VECGEOM_NAVIGATOR
#define USE_VECGEOM_NAVIGATOR
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SimplePhysicsList::SimplePhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SimplePhysicsList::~SimplePhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimplePhysicsList::ConstructParticle()
{

// G4EmStandardPhysics particles:: 
  // gamma
  G4Gamma::Gamma();

  // leptons
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();

  // mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();

  // barions
  G4Proton::Proton();
  G4AntiProton::AntiProton();

  // ions
  G4Deuteron::Deuteron();
  G4Triton::Triton();
  G4He3::He3();
  G4Alpha::Alpha();
  G4GenericIon::GenericIonDefinition();

// G4DecayPhysics particles::
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  

// G4HadronPhysicsFTFP_BERT and G4HadronElasticPhysics particles::
//  G4MesonConstructor pMesonConstructor;
//  pMesonConstructor.ConstructParticle();

//  G4BaryonConstructor pBaryonConstructor;
//  pBaryonConstructor.ConstructParticle();

//  G4ShortLivedConstructor pShortLivedConstructor;
//  pShortLivedConstructor.ConstructParticle();  

//  G4IonConstructor pConstructor;
//  pConstructor.ConstructParticle();  


  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 
/*
std::cout<< "BEGIN THE =" <<G4ParticleTable::GetParticleTable()->size() <<std::endl;   
  TPartIndex      *pIndex       = TPartIndex::I();
  G4ParticleTable *thePartTable = G4ParticleTable::GetParticleTable();
  
  G4int tabNumPart = pIndex->NPart();
  for(G4int ip = 0; ip < tabNumPart; ++ip){
     G4int pdg = pIndex->PDG(ip);
     G4ParticleDefinition *partDef = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
     if(!partDef)
     
  }
*/
std::cout<< "END THE =" <<G4ParticleTable::GetParticleTable()->size() <<std::endl;   
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimplePhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructTotal();
  // ConstructDecay();  // Should be embedded in the 'Total' Process - tbc

  //  HadronPhysicsFTFP_BERT_WP();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimplePhysicsList::ConstructTotal()
{
  // G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  G4cout << " SimplePhysicsList: constructing one TotalPhysicsProcess per particle " << G4endl;
  
  // Must use the old functionality for adding processes
  //  - the new one works only for recognised processes, e.g. Brem, compton ..
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    // G4String particleName = particle->GetParticleName();
    
    G4VRestContinuousDiscreteProcess* totalPhysics=
       new TotalPhysicsProcess("TabulatedPhysics");

    // ph->RegisterProcess( totalPhysics, particle);
    pmanager->AddContinuousProcess( totalPhysics );
    pmanager->AddDiscreteProcess( totalPhysics );
    pmanager->AddRestProcess( totalPhysics );

    // Add UserSpecialCut
    //pmanager->AddProcess(new G4UserSpecialCuts, 1, 1, 1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
#include "G4Decay.hh"

void SimplePhysicsList::ConstructDecay()
{
  // Add Decay Process
  G4cout << "Constructing separate Decay process(es) for each particle." << G4endl;
  
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    
    if (theDecayProcess->IsApplicable(*particle))
    {
      pmanager ->AddProcess(theDecayProcess);
      
      // set ordering for PostStepDoIt and AtRestDoIt
      
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}


*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void SimplePhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "SimplePhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue // / CLHEP::Unit::mm
    << " mm" << G4endl;
    // G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "proton");

  if (verboseLevel>0) DumpCutValuesTable();
}
*/
/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimplePhysicsList::HadronPhysicsFTFP_BERT_WP()
{
}
*/
