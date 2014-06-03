#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIExecutive.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "CLHEP/Random/RanluxEngine.h" 
#include "MyDetectorConstruction.hh" 
#include "MyPrimaryGeneratorAction.hh" 
#include "MyEventAction.hh" 
#include "MyUserActionInitialization.hh"
// --- Physics Lists ---
#include "G4PhysListFactory.hh"
// --- 

//01.25.2009 Xin Dong: This example came from the original sequential
//program FullCMS. The original program is changed here to support parallel
//computing with multiple threads. All events are assigned to each worker
//thread in a round robin fashion. All threads share most detector data 
//including physics table and physics vector for some physics processes.
//The master process initializes the data in a regular way. However, worker
//threads initialize thread private data only.
//#include "G4MTParTopC.icc"

//01.25.2009 Xin Dong: Threads share this object.
//MyDetectorConstruction *detector = 0;

int main(int argc,char** argv) { 
  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  G4Random::setTheEngine( &defaultEngine );
  G4int seed = 1220515164;
  G4Random::setTheSeed( seed );

  G4cout << G4endl 
         << " ===================================================== " << G4endl 
         << " Initial seed = " << seed << G4endl 
	 << " ===================================================== " << G4endl 
	 << G4endl; 
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager; // = new G4RunManager; 
  //G4int nt = 2;
  // GF: note argv[1] : input file in batch
  //              [2] : physics lists to use - optional
  //if ( argc > 1 ) nt = atoi( argv[1] ); 
  runManager->SetNumberOfThreads( 2 );
#else
  G4RunManager* runManager = new G4RunManager;
#endif
  //01.25.2009 Xin Dong: The master thread and worker threads have different behavior
  //in the phase of initialization.
  //if (threadRank == 0)
  //  runManager = new G4RunManager;
  //else
  //  runManager = new G4RunManager(1);

  //01.25.2009 Xin Dong: The master thread and worker threads have different behavior.
  //in the phase of initialization.
  //if (threadRank == 0)
    MyDetectorConstruction* detector = new MyDetectorConstruction;
    //else
    //detector->SlaveMyDetectorConstruction();

  runManager->SetUserInitialization(detector); 

  // --- Physics Lists ---
  G4String namePL;
  if ( argc > 2 ) { // The second argument, when present, is a Physics List.
    namePL = argv[2];
    G4cout << "ParFullCMS test using physics list " << namePL << G4endl;
  } else {
    namePL = "QGSP_BERT";
    G4cout << "ParFullCMS using physics list QGSP_BERT (default)" << G4endl;
  }

  G4PhysListFactory factory;
  if ( factory.IsReferencePhysList( namePL ) ) {
    runManager->SetUserInitialization( factory.GetReferencePhysList( namePL ) );
  } else {
    G4cerr << "ERROR: Physics List " << namePL << " UNKNOWN!" << G4endl;
  }

  //thePL->SetDefaultCutValue( 0.020 *mm ); // 20 microns 

  runManager->SetUserInitialization( new MyUserActionInitialization );


#ifdef G4VIS_USE
   // Visualization manager
   //
   G4VisManager* visManager = new G4VisExecutive;
   visManager->Initialize();
#endif

  runManager->Initialize(); 

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
                                                                                
  if ( argc==1 ) {   // Define UI session for interactive mode.
      G4UIExecutive* ui = new G4UIExecutive(argc,argv);
    //G4UIsession * session = new G4UIterminal(new G4UItcsh);
    UI->ApplyCommand("/control/execute vis.mac");
    G4cout << "Now, please, apply beamOn command..." << G4endl;
      ui->SessionStart();
    //session->SessionStart();
    delete ui;
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

#ifdef G4VIS_USE
  delete visManager;
#endif

  //01.25.2009 Xin Dong: The master thread and worker threads have different behavior
  //in the phase of exit.
  //if (threadRank != 0) detector->SlaveDestroy();

  //01.25.2009 Xin Dong: The master thread and worker threads have different behavior
  //in the phase of exit. It is better that all threads destroy this object. However,
  //it results in panics by destroying shared data many times.
  //if (threadRank == 0) 
  //  delete runManager;
  return 0; 
} 
