//--------------------------------------------------------
// Last update: 18-Mar-2010
//
// The Physics List can be specified in the command line,
// as the second argument, e.g. :
//
//          mainHadCalo  run1.g4  FTF_BIC
//
// Notice that the name of the Physics List must be in
// upper case, exactly as the corresponding header file.
// If the Physics List is not specified, then  QGSP_BERT
// will be used by default.
// 
//--------------------------------------------------------

#include "StatAccepTestDetectorConstruction.hh" 
#include "StatAccepTestPrimaryGeneratorAction.hh" 
#include "StatAccepTestEventAction.hh" 
#include "StatAccepTestRunAction.hh" 
#include "StatAccepTestTrackingAction.hh" 
#include "StatAccepTestAnalysis.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "G4RunManager.hh" 
#include "G4UImanager.hh" 
#include "G4UIterminal.hh" 
#ifdef G4UI_USE_XM 
#include "G4UIXm.hh" 
#endif 
#include "CLHEP/Random/RanluxEngine.h" 
#include <ctime> 


int main( int argc, char** argv ) { 

  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  G4int seed = time( NULL ); 
  CLHEP::HepRandom::setTheSeed( seed ); 
  G4cout << G4endl 
         << " ===================================================== " << G4endl 
         << " Initial seed = " << seed << G4endl 
	 << " ===================================================== " << G4endl 
	 << G4endl; 
  G4RunManager* runManager = new G4RunManager; 
  
  G4String namePL;
  if ( argc > 2 ) { // The second argument, when present, is a Physics List.
    namePL = argv[2];
    G4cout << "HadCalCMS test using physics list " << namePL << G4endl;
  } else {
    namePL = "QGSP_BERT";
    G4cout << "HadCalCMS using physics list QGSP_BERT (default)" << G4endl;
  }

  G4PhysListFactory factory;
  if ( factory.IsReferencePhysList( namePL ) ) {
    G4VModularPhysicsList* physList = factory.GetReferencePhysList( namePL );
    runManager->SetUserInitialization( physList );
  } else {
    G4cerr << "ERROR: Physics List " << namePL << " UNKNOWN!" << G4endl;
  }

  runManager->SetUserInitialization( new StatAccepTestDetectorConstruction ); 
  runManager->SetUserAction( new StatAccepTestPrimaryGeneratorAction ); 
  runManager->SetUserAction( new StatAccepTestRunAction ); 
  runManager->SetUserAction( new StatAccepTestEventAction ); 
  runManager->SetUserAction( new StatAccepTestTrackingAction ); 
  runManager->Initialize(); 

  G4UImanager* UI = G4UImanager::GetUIpointer(); 

  if ( argc == 1 ) {   // Define UI session for interactive mode. 
    G4UIsession* session = 0; 
#ifdef G4UI_USE_XM 
    session = new G4UIXm(argc,argv); 
#else 
    session = new G4UIterminal(); 
#endif 
#ifdef G4UI_USE_XM 
    // Customize the G4UIXm menubar with a macro file : 
    UI->ApplyCommand( "/control/execute gui.g4" ); 
#else 
    G4cout << "Now, please, apply beamOn command..." << G4endl; 
#endif 
    session->SessionStart(); 
    delete session; 
  } else {   // Batch mode 
    G4String command = "/control/execute "; 
    G4String fileName = argv[1]; 
    UI->ApplyCommand( command + fileName ); 
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
