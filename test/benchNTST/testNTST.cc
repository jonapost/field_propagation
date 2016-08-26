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
// $Id: testNTST.cc,v 1.3 2007-10-26 09:51:28 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------
//      GEANT 4 - exampleNTST : BaBar SVT standalone simulation
//
//      For information related to this code contact:
//         Tatiana Nikitina - Tatiana.Nikitina@cern.ch
// ----------------------------------------------------------------

#include <sstream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4String.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "NTSTFileRead.hh"
#include "NTSTDetectorConstruction.hh"
#include "NTSTPrimaryGeneratorAction.hh"
#include "NTSTRunAction.hh"
#include "NTSTEventAction.hh"
#include "NTSTSteppingAction.hh"

// --- Physics Lists ---

#include "NTSTPhysicsList.hh"

#if 0 // def flagLHEP
#include "LHEP.hh"
#endif

#if 0 // def flagQGSP_EMV
#include "QGSP_EMV.hh"
#endif

// #ifdef flagFTFP
#if 1
#include "FTFP_BERT.hh"
#endif

#if 1
// #ifdef flagQGSP_BERT
#include "QGSP_BERT.hh"
#endif

// #if 0 // def flagQBBC
#if 1
#include "QBBC.hh" 
#endif

// #ifdef flagQGSP_BIC 
#if 1
#include "QGSP_BIC.hh" 
#endif

#if 0 // def flagQGSP_BERT_HP
#include "QGSP_BERT_HP.hh"
#endif

//  The new way of configuring different physics lists - copied from Hadr00
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

// --- 

#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"

int main(int argc,char** argv)
{

   // Construct the default run manager
   G4RunManager * runManager = new G4RunManager;
    
   // set mandatory initialization classes
   NTSTDetectorConstruction* detector = new NTSTDetectorConstruction();
     
   // detector construction object
   runManager->SetUserInitialization(detector);

   // physics list

#ifdef OLD_PHYSICS_LIST
   G4VUserPhysicsList* thePL = 0;
#ifdef flagQGSP
  thePL = new QGSP;   
#else
  #ifdef flagQGSP_EMV
    thePL = new QGSP_EMV;
  #else
    #ifdef flagQBBC
      thePL = new QBBC;
    #else
      #ifdef flagFTFP
        thePL = new FTFP_BERT;
      #else
        #ifdef flagQGSP_BERT
          thePL = new QGSP_BERT;
        #else
          #ifdef flagQGSP_BIC
            thePL = new QGSP_BIC;
          #else
            thePL = new NTSTPhysicsList();  // Default !
          #endif
        #endif
      #endif
    #endif
  #endif
#endif
  runManager->SetUserInitialization( thePL ); 

#else
  //  The new way of configuring different physics lists - copied from Hadr00
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = 0;
  G4String physName = "";

  // Physics List name defined via 3nd argument
  if (argc>=3) { 
    physName = argv[2];
    G4cout << "Physics List from argv[2]:  requested " 
           << physName << G4endl;
  }

  // Physics List is defined via environment variable PHYSLIST
  if("" == physName) {
    char* path = getenv("PHYSLIST");
    if (path) { 
       physName = G4String(path);
       G4cout << "Physics List from PHYSLIST environment variable: requesting " 
           << physName << G4endl; 
    }
  }

  // if name is not known to the factory use FTFP_BERT
  if("" == physName || !factory.IsReferencePhysList(physName)) {
     G4cout << "Warning> Physics List choice '" << physName <<  "' was not found. " << G4endl; 
     physName = "FTFP_BERT"; 
     G4cout << "Using default   Physics List " << physName << G4endl;
  }else{
     G4cout << "Using requested Physics List " << physName << G4endl;
  }

  // reference PhysicsList via its name
  phys = factory.GetReferencePhysList(physName);
  if ( phys == 0 ) { 
     G4cerr << "ERROR> Request to factory to create Physics List '" << physName <<  "' failed. " ; 
     exit(1); 
  } 
  runManager->SetUserInitialization(phys);
#endif

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new NTSTPrimaryGeneratorAction);
  runManager->SetUserAction(new NTSTRunAction);
  runManager->SetUserAction(new NTSTEventAction);
  runManager->SetUserAction(new NTSTSteppingAction);
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI terminal for interactive mode  
    { 
      G4UIsession * session = new G4UIterminal(new G4UItcsh);
      //      UI->ApplyCommand("/control/execute prerunNTST.mac");    
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

