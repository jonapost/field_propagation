//--------------------------------------------------------
// Last update: 12-Nov-2009
//
// The Physics List can be specified in the command line,
// as the second argument, e.g. :
//
//          mainApplication  bench1.g4  FTF_BIC
//
// Notice that the name of the Physics List must be in
// upper case, exactly as the corresponding header file.
// If the Physics List is not specified, then  QGSP_BERT
// will be used by default.
//
//--------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "Randomize.hh" 
#include "MyDetectorConstruction.hh" 
#include "MyPrimaryGeneratorAction.hh" 
#include "MyEventAction.hh" 
#include <ctime> 


int main(int argc,char** argv) { 

  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  G4Random::setTheEngine( &defaultEngine ); 
  G4int seed = time( NULL ); 
  G4Random::setTheSeed( seed ); 
  G4cout << G4endl 
         << " ===================================================== " << G4endl 
         << " Initial seed = " << seed << G4endl 
	 << " ===================================================== " << G4endl 
	 << G4endl; 

  G4RunManager* runManager = new G4RunManager; 

  G4String namePL;
  if ( argc > 2 ) { // The second argument, when present, is a Physics List.
    namePL = argv[2];
    G4cout << "FullCMS test using physics list " << namePL << G4endl;
  } else {
    namePL = "QGSP_BERT";
    G4cout << "FullCMS using physics list QGSP_BERT (default)" << G4endl;
  }

  G4PhysListFactory factory;
  if ( factory.IsReferencePhysList( namePL ) ) {
    G4VModularPhysicsList* physList = factory.GetReferencePhysList( namePL );
    runManager->SetUserInitialization( physList );
  } else {
    G4cerr << "ERROR: Physics List " << namePL << " UNKNOWN!" << G4endl;
  }

  runManager->SetUserInitialization( new MyDetectorConstruction ); 
  runManager->SetUserAction( new MyPrimaryGeneratorAction ); 
  runManager->SetUserAction( new MyEventAction ); 

  runManager->Initialize(); 

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
                                                                                
  if ( argc==1 ) {   // Define UI session for interactive mode.
    G4UIsession * session = new G4UIterminal();
    UI->ApplyCommand("/control/execute vis.mac");
    G4cout << "Now, please, apply beamOn command..." << G4endl;
    session->SessionStart();
    delete session;
  } else {   // Batch mode 
    G4String command = "/control/execute "; 
    G4String fileName = argv[1]; 
    UI->ApplyCommand(command+fileName); 
  } 

  G4cout << G4endl 
	 << " ===================================================== " << G4endl 
         << " Final random number = " 
         << CLHEP::HepRandom::getTheEngine()->flat() << G4endl 
	 << " ===================================================== " << G4endl 
         << G4endl; 

  delete runManager;

  return 0; 
} 
