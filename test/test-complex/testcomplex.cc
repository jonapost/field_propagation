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
//  Complex geometry test as a modified Geant4 example.
//  Original Geant4 example: examples/extended/persistency/gdml/G01
//
#include <vector>

#include "G4RunManager.hh"
#include "G4UImanager.hh"

// incs. for adding a region for each logical volume
#include "G4LogicalVolumeStore.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"

// incs. for setting the production cuts
#include "G4RToEConvForGamma.hh"
#include "G4RToEConvForElectron.hh"
#include "G4RToEConvForPositron.hh"
#include "G4RToEConvForProton.hh"
#include "G4MaterialCutsCouple.hh"

// incs. for setting up tracking cuts for gamma,e-,e+,p+
#include "G4UserLimits.hh"

// for adding UserSpecialCut process 
#include "G4UserSpecialCuts.hh"

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"

#if GEANT4_LESS_10 == 1
#include "G4StepLimiterBuilder.hh"
#else
#include "G4StepLimiterPhysics.hh"
#endif

#include "QBBC.hh"
#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4GDMLParser.hh"


#include "SimplePhysicsList.hh"
#include "TabulatedDataManager.hh"
#include "TotalPhysicsProcess.hh"


#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

#include <getopt.h>

#include "CMSApp.hh"
#include "LogicalVolumeHandler.hh"

///////////////////////////////////////////////////////////////////////////////////////
void AddRegionsToLogicalVols(G4double trackingCutInEnergy);
void SetG4ProductionCut(G4double ethresh);
void usage();
void parsearg(int argc, char **argv, std::string &geomFile, std::string &eventFile,
              std::string &g4macroFile, G4double &cutvalue, std::string &physListName,
              int &scoreType, bool &useUniformField);
////////////////////////////////////////////////////////////////////////////////////////


int main(int argc,char **argv)
{

  std::string geomFile     = "";        // GDML geometry file
  std::string eventFile    = "";        // root file with pre-generated primary events
  std::string g4macroFile  = "";        // Geant4 macro file
  std::string physListName = "TABPHYS"; // physics list name
  G4double cutvalue        = 0.;        // common tracking and production cut [GeV]
  G4int scoreType          = 0;         // score type flag
  bool useUniformField = false;

  //
  // Parsing arguments
  //
  parsearg(argc, argv,  geomFile, eventFile, g4macroFile, cutvalue, physListName,
           scoreType, useUniformField);
  //
  // Print information 
  //
  G4cout<< G4endl <<
  "-----------------------------------------------------------------------------"
  << G4endl <<
  "|  testcomplex was executed:                                                 "
  << G4endl <<
  "|         Geometry file       : "<< geomFile               << G4endl <<
  "|         Magnetic field      : "<< ( useUniformField ? "Uniform" : "Interpolated") << G4endl <<
//  "|         Magnetic field file : "<< magFieldFile           << G4endl <<
  "|         Primary events      : "<< eventFile              << G4endl <<
  "|         Geant4 macro        : "<< g4macroFile            << G4endl <<
  "|         Physics list        : "<< physListName           << G4endl <<
  "|         Low energy cut[GeV] : "<< cutvalue               << G4endl <<
  "|         Score type          : "<< scoreType              << G4endl <<
  "|         Xsec data file      : "<< getenv("VP_DATA_XSEC") << G4endl <<
  "|         Fstate data file    : "<< getenv("VP_DATA_FSTA") << G4endl <<
  "-----------------------------------------------------------------------------"
  << G4endl << G4endl;

  bool useUI = false;
  // bool useUI= (g4macrofile == std::string("-");

  // bool useParticleGun= (eventFile == std::string("-");
  
  G4double tabPhysEnergyLimit     = cutvalue;              // low energy cut in tab. physics [GeV]
  G4double fTrackingCutInEnergy   = cutvalue * CLHEP::GeV; // same used in Geant4
  G4double fProductionCutInEnergy = cutvalue * CLHEP::GeV; // same used in Geant4 as production cut

  //
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);
  // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);  

  //
  // Get the geometry from the GDML file
  //
  G4GDMLParser parser;
  //parser.SetStripFlag(FALSE);
  parser.Read(geomFile);

  //
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;


  //
  // Set mandatory initialization classes:
  //   DetectorConstruction, PhysicsList
  //
  DetectorConstruction* dc=
     new DetectorConstruction( parser.GetWorldVolume(), useUniformField );

  runManager->SetUserInitialization( dc );

  // WITH G4 10.x
  if(std::string("TABPHYS") == physListName) {
    runManager->SetUserInitialization(new SimplePhysicsList);
    printf("Physics: Using TABPHYS - ie SimplePhysicsList");
    RunAction::isTabPhys = TRUE;
  } else {
    RunAction::isTabPhys = FALSE;
    G4VModularPhysicsList* physicsList = 0;
    if(std::string("QBBC") == physListName) {
      physicsList = new QBBC;
      printf("Physics: Using      QBBC    Physics List");      
    } else if(std::string("FTFP_BERT") == physListName) {
      physicsList = new FTFP_BERT;
      printf("Physics: Using FTFP_BERT    Physics List");            
    } else if(std::string("FTFP_BERT_HP") == physListName) {
      physicsList = new FTFP_BERT_HP;
      printf("Physics: Using FTFP_BERT_HP Physics List");                  
    } else {
      G4cout << G4endl <<
      "**************************************************************************"
      << G4endl <<  
      "*                          Unknown physics list                          *"           
      << G4endl <<  
      "**************************************************************************"
      << G4endl << G4endl;  
      usage(); 
      exit(EXIT_FAILURE);
    }
#if GEANT4_LESS_10 == 1
      physicsList->RegisterPhysics(new G4StepLimiterBuilder());
#else 
      physicsList->RegisterPhysics(new G4StepLimiterPhysics());
#endif
      runManager->SetUserInitialization(physicsList); 
  }

  //
  // Setting some variables that control the tracking 
  // 
    //TabulatedDataManager::SetVerboseLevel( 4 );
  TotalPhysicsProcess::SetEnergyLimit( tabPhysEnergyLimit ); 
  TabulatedDataManager::fgIsUseRange = TRUE;
  SteppingAction::fgTrackingCutInEnergy = fTrackingCutInEnergy;  
  //G4SteppingManager::fgProductionCutInEnergy = fProductionCutInEnergy;

  CMSApp::fgIsScoreActive = FALSE; 
  if(scoreType > 1)
    CMSApp::fgIsScoreActive = TRUE;
 
  RunAction::fgScoreTypeFlag = scoreType;
  // 
  // Set user action classes
  //
  runManager->SetUserAction(new PrimaryGeneratorAction(eventFile));
  //
  runManager->SetUserAction(new RunAction);
  //
  runManager->SetUserAction(new EventAction);
  //
  runManager->SetUserAction(new SteppingAction);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4String command_ls = "ls /gun";
  UImanager->ApplyCommand(command_ls);

  G4String command_query = "cd /gun ; ?particle";
  UImanager->ApplyCommand(command_query);
  
  //
  // Initialize Geant4 run manager
  // 
  runManager->Initialize();
  G4cout<< "----> Run manager has been initialized. " <<G4endl;

  // 
  // Check logical volumes that are used in teh geometry 
  //
  LogicalVolumeHandler::Instance(); // just for initialization

 
  // 
  // When running with Geant4 physics list:
  //   -loop over the LogicalVolumes and define a G4Region for each of them
  //   -add a tracking cut in energy as UserSpecialCut to each region  
  //   -set production cuts in energy at each region defined above 
  //
  if(std::string("TABPHYS") != physListName) {
    AddRegionsToLogicalVols(fTrackingCutInEnergy);
    SetG4ProductionCut(fProductionCutInEnergy);
  }

  //
  // Execute the given Geant4 macro files
  // 
  if(std::string("TABPHYS") != physListName) {
    G4cout << G4endl <<
    "----> Initialiation of Geant4 physics tables can take several minutes in case of high \n"
    "      number of regions. Physics tables (dEdx, lambda) are built for each region but  \n"
    "      progress is reported below only for one of them (i.e.DefaultRegionForTheWorld)  \n"
    "      ......You have time for a coffee....!" << 
    G4endl;  
  }
  // G4UImanager* UImanager = G4UImanager::GetUIpointer();

  G4UIExecutive* ui = 0;
  if ( useUI ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4cout << "Batch mode - no UI exec" << G4endl;
    G4String command = "/control/execute ";
    G4cout << "Macro file name= " << g4macroFile << G4endl;
     
    UImanager->ApplyCommand(command+g4macroFile);
  }  
  else {
    // interactive mode
    ui->SessionStart();
    delete ui;
  }
  
  //
  // End of story
  //
  delete runManager;

  return 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void AddRegionsToLogicalVols(G4double trackingCutInEnergy){
   G4ProductionCuts pcut;
   pcut.SetProductionCut(1.0*CLHEP::mm); // will be set properly later 

   G4cout<< G4endl<<
   "------------------- Adding regions to logical volumes ---------------------"
   << G4endl;

   G4Region *region = 0;
   const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
   int  numlvs = lvs->size();  // number of logical volumes 
   G4cout<< "---- Number of logical volumes =  "<< lvs->size()<< G4endl;
   const G4RegionStore* rgs  =  G4RegionStore::GetInstance();
   G4cout<< "---- Number of current regions =  "<< rgs->size() <<G4endl;
   std::vector<G4Region*>::const_iterator rgciter;
   for( rgciter = rgs->begin(); rgciter != rgs->end(); rgciter++ )
      G4cout<<"        --- region name = "<< (*rgciter)->GetName()<< G4endl;
  
   //
   // Handle root logical volume that belongs to DefaultRegionForTheWorld   
   //
   //G4cout<< "# volumes in DefaultRegionForTheWorld  = "<<(*(rgs->begin()))->GetNumberOfRootVolumes()<<G4endl;
   std::vector<G4LogicalVolume*>::const_iterator lviter= (*(rgs->begin()))->GetRootLogicalVolumeIterator();
   G4LogicalVolume *lgvWorld = *lviter;
   G4Region *regWorld = *(rgs->begin());   

   char regname[14];
   //
   // Loop over logical volumes and define a region for each of them and add a UserLimit 
   // i.e. a tracking cut in energy. Don't add region to the root logical volume of world, becuse it already
   // has a DefaultRegionForTheWorld. 
   //
   for( int lvindx = 0; lvindx < numlvs; ++lvindx )
   {
     G4LogicalVolume *glv = (*lvs)[lvindx];
//     if(lvindx==4683) continue;   // THIS is DefaultRegionForTheWorld
     if(glv == lgvWorld)
     {
       G4cout << "The "<< lvindx << " logical volume is the world that has region by default."<<G4endl;  
       regWorld->SetProductionCuts(new G4ProductionCuts(pcut));
       regWorld->RegisterMaterialCouplePair( glv->GetMaterial(),
                                         new G4MaterialCutsCouple( glv->GetMaterial(), new G4ProductionCuts(pcut))
                                       );
       // setting user limits to be able to set tracking cuts
       glv->SetUserLimits( new G4UserLimits( DBL_MAX, DBL_MAX, DBL_MAX, trackingCutInEnergy, 0.0 ) );
       //regWorld->SetUserLimits( new G4UserLimits( DBL_MAX, DBL_MAX, DBL_MAX, trackingCutInEnergy, 0.0 ) );
       continue;
     }

     snprintf(regname,13,"Region_%4.4d_",lvindx);
     region = new G4Region(G4String(regname)+glv->GetName());
     region->AddRootLogicalVolume(glv);
     region->SetProductionCuts(new G4ProductionCuts(pcut));
     region->RegisterMaterialCouplePair( glv->GetMaterial(),
                                         new G4MaterialCutsCouple( glv->GetMaterial(), new G4ProductionCuts(pcut))
                                       );
     // setting user limits to be able to set tracking cuts
     glv->SetUserLimits( new G4UserLimits( DBL_MAX, DBL_MAX, DBL_MAX, trackingCutInEnergy, 0.0 ) );
     //region->SetUserLimits( new G4UserLimits( DBL_MAX, DBL_MAX, DBL_MAX, trackingCutInEnergy, 0.0 ) );
   }
   G4cout<< "---- Number of current regions at the end =  "<< rgs->size() <<G4endl;
   G4cout<< "---- User limits were also set for each logical volume in order to \n"
            "     be able to use tracking cut (in energy)."<< G4endl;
   G4cout<<
   "---------------------------------------------------------------------------"
   << G4endl;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// Setting production cut (for gamma, e-, e+, p+) to a value given in energy.
// Will be done only if one of the original G4 physics list is used.
void SetG4ProductionCut(G4double ethresh){
  //G4double ethresh = 0.001*CLHEP::GeV; // production cut in energy to be set in range 1 MeV now
   
  G4cout <<
  "============================================================================"
  << G4endl;
  printf("\t\tSetting Geant4 production cuts in energy START\n");
  G4MaterialTable *theMaterialTable = (G4MaterialTable*)G4Material::GetMaterialTable();
  //printf("Number of materials = %d\n",theMaterialTable->size());

  G4VRangeToEnergyConverter *converter[4];
  converter[0] = new G4RToEConvForGamma();
  converter[1] = new G4RToEConvForElectron();
  converter[2] = new G4RToEConvForPositron();
  converter[3] = new G4RToEConvForProton();
   
  
   
  G4RegionStore *theRegionStore = (G4RegionStore*) G4RegionStore::GetInstance();
      
  G4int nRegions = theRegionStore->size();
  printf("\t\tWe found %d regions:\n",nRegions);

  for(G4int i=0;i<nRegions;++i) {
    G4Region *reg = (*theRegionStore)[i];
//    std::cout<<"\t\tName of "<< i <<"-th region is: "<< reg->GetName() << std::endl;

    for(unsigned int imat=0; imat<theMaterialTable->size(); ++imat) {
      G4MaterialCutsCouple *mcCouple = reg->FindCouple((*theMaterialTable)[imat]);
      if(mcCouple!=0) {
//         G4ProductionCuts* pcuts = reg->GetProductionCuts();
//         const char* parcuts[4]={"gamma","e-","e+","proton"};
//         printf("Production Cuts (cm) : ");
//         for(G4int ic=0; ic<4; ++ic) 
//           printf("%s=%12.2g  [cm]   ",parcuts[ic],pcuts->GetProductionCut(ic)/CLHEP::cm);
//         printf("\n");
         G4Material *mat = (*theMaterialTable)[imat];  
//         printf("Production Cuts (GeV): ");
//         for(G4int ic=0; ic<4; ++ic)
//           printf("%s=%12.2g  [GeV]  ", parcuts[ic],
//                converter[ic]->Convert(pcuts->GetProductionCuts()[ic],mat)/CLHEP::GeV);
//         printf("\n");

         for(G4int ic=0; ic<4; ++ic) {
           G4double lmin = 1.0*CLHEP::nm;
           G4double lmax = 10.0*CLHEP::km;
           while(std::abs(lmin-lmax) > 0.5*(lmin+lmax)*1e-10) {
             G4double lmid = 0.5*(lmin+lmax);
             if((converter[ic]->Convert(lmid,mat))/(CLHEP::GeV) > ethresh/(CLHEP::GeV)) lmax=lmid;
             else lmin = lmid;
           }
           reg->GetProductionCuts()->SetProductionCut(0.5*(lmin+lmax),ic);
         }

         //After setting
//         printf("-------------------: After Setting them ----------------------\n");
//         pcuts = reg->GetProductionCuts();
//         printf("Production Cuts (cm) : ");
//         for(G4int ic=0; ic<4; ++ic) 
//           printf("%s=%12.2g  [cm]   ",parcuts[ic],pcuts->GetProductionCut(ic)/CLHEP::cm);
//         printf("\n");
//         printf("Production Cuts (GeV): ");
//         for(G4int ic=0; ic<4; ++ic)
//           printf("%s=%12.2g  [GeV]  ", parcuts[ic],
//         converter[ic]->Convert(pcuts->GetProductionCuts()[ic],mat)/CLHEP::GeV);
//         printf("\n#############################################################################\n");
      }
    }
  }

  printf("\t\tSetting Geant4 production cuts in energy END\n");
  G4cout <<
  "============================================================================"
  << G4endl;

  // delete converters     
  for(G4int i=0; i<4; ++i) delete converter[i]; 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void usage()
{
  G4cerr <<  G4endl <<
  "============================================================================"
  << G4endl <<
  "[NAME]" << G4endl <<
  "    complextest ----------- Geant4 test with complex geometry ------------- " 
       << G4endl << G4endl <<
  "[SYNOPSIS]" << G4endl <<
  "    testcomplex --geomFile [FILE] --eventFile [FILE] --g4macroFile [FILE] \n"
  "                --lowEnergyCut [VALUE] < --physListName [NAME] >          \n"
  "                < --uniformField | --bilinearField  >                       "
        << G4endl << G4endl <<
  "[DESCRIPTION]" << G4endl <<
  "    Run Geant4 on a complex geometry setup given in form of GDML input file   \n"
  "    (cms2015.gdml) using some pre-generated events given in form of root file \n"
  "    (14TeVminbias.root) with optional physics list. The default physics list  \n"
  "    is the tabulated physics that needs the cross sections and final states   \n"
  "    files given by setting the VP_DATA_XSEC=where/your/xsec_*.root and        \n"
  "    VP_DATA_FSTA=where/your/fstate_*.root environment variables. Note, that   \n"
  "    these data correspond to a given cut value in energy that you need to     \n"
  "    provide to testcomplex as an input parameter through the -l flag in [GeV] \n"
  "    units (this value is indicated at the end of the xsec/fstate_*.root files). "
        << G4endl << G4endl <<
  "[REQUIRED]" << G4endl <<
  " --geomFile     [FILE] :  geometry GDML file i.e. where/your/cms2015.gdml file\n"
  "                          is located                                            " 
          << G4endl <<  
  " --eventFile    [FILE] : primary event generator file i.e. where/your/pp14TeVm\n"
  "                         inbias.root file is located.                           "
          << G4endl <<
  " --geantMacro   [FILE] : Geant4 macro file i.e. g4mac.mac                       " 
          << G4endl <<
  " --lowEnergyCut [VALUE]: low energy cut value in the tabulated physics [GeV]    " 
        << G4endl << G4endl <<

  "[OPTIONAL]" << G4endl <<
  " --physListName <NAME> : physics list: TABPHYS (default), FTFP_BERT, FTFP_BERT\n"
  "                         _HP, QBBC                                              "
          << G4endl <<
  " --uniformField        : use uniform magnetic field                             "
          << G4endl <<
  " --bilinearField       : use bilinear CMS magnetic field                        "
          << G4endl <<
  // " --useUniform   [FLAG] : whether to use uniform magnetic field                  "
  //       << G4endl <<
  // "                         Potential values: y/Y/t/T/1 = true, n/N/f/F/0 = false  "
  //      << G4endl <<
  " --scoreType    [VALUE]: type of scoring:                                     \n"
  "                          0 : no scoring, only run-time is reported (default) \n"
  "                          1 : option 1 plus step statistics                   \n"
  "                          2 : option 2 plus histograms in CMS ECAL            \n"
          << G4endl <<
  " --info                : print all information regarding testcomplex          \n"
  << G4endl <<
  "============================================================================"
  << G4endl << G4endl;
}

void parsearg(int argc, char **argv, std::string &geomFile, std::string &eventFile,
              std::string &g4macroFile, G4double &cutvalue, std::string &physListName,
              int &scoreType, bool &useUniformField){
    if(argc != 2 && argc < 5) {
      G4cout << G4endl <<
      "**************************************************************************"
      << G4endl <<  
      "*                        Missing required parameters                     *"           
      << G4endl <<  
      "**************************************************************************"
      << G4endl << G4endl;  
      usage(); 
      exit(EXIT_FAILURE);
    }

    int opt= 0;
    int counter   = 0;
     
    static struct option long_options[] = {
        {"geomFile"      ,  required_argument, 0,  0 },
        {"eventFile"     ,  required_argument, 0,  0 },
        {"geantMacro"    ,  required_argument, 0,  0 },
        {"lowEnergyCut"  ,  required_argument, 0,  0 },
        {"physListName"  ,  required_argument, 0,  0 },
        {"scoreType"     ,  required_argument, 0,  0 },
        {"uniformField"  ,  no_argument      , 0,  0 },
        {"bilinearField" ,  no_argument      , 0,  0 },
        {"info"          ,  no_argument      , 0,  0 },
        {0               ,  0                , 0,  0 }
    };

    int long_index =0;
    while ((opt = getopt_long_only(argc, argv,"", long_options, &long_index )) != -1) {
        switch (opt) {
             case  0  : 
                        switch (long_index) {
                          case  0  : geomFile     = optarg; 
                                     ++counter;
                              break;
                          case  1  : eventFile    = optarg;
                                     ++counter;
                              break;
                          case  2  : g4macroFile  = optarg;
                                     ++counter;
                              break;
                          case  3  : cutvalue     = atof(optarg);  
                                     ++counter;
                              break;
                          case  4  : physListName = optarg; 
                              break;
                          case  5  : scoreType    = atof(optarg);
                              break;
                          case  6  : useUniformField = true;
                                     G4cout << " Chose Uniform field " << G4endl;                             
                              break;
                          case  7  : useUniformField = false;   // bilinear field
                                     G4cout << " Chose Bilinear field " << G4endl;
                              break;                              
                              /****
                                     if( optarg[0] == 'y' || optarg[0] == 'Y'
                                         || optarg[0] == 't' || optarg[0] == 'T' )
                                     {
                                       useUniformField = true;
                                     } else {
                                        if( optarg[0] == 'n' || optarg[0] == 'N'
                                            || optarg[0] == 'f' || optarg[0] == 'F' )
                                        {
                                           useUniformField = true;
                                        } else {
                                           if( optarg[0] >= '0' && optarg[0] <= '9' ) {
                                              useUniformField = ( optarg[0] > '0' );
                                           } else {
                                              std::cerr << " Invalid value: useUniform = " << optarg
                                                        << "  Expected: y/Y/t/T/1 = true, n/N/f/F/0 = false"; 
                                              useUniformField = true;
                                           }
                                        }
                                     }
                              break;
                              ****/
                          case  8  : usage();
                                     exit(EXIT_SUCCESS); 
                              break;
                        }
                 break; 
             case '?' : usage();
                 exit(EXIT_FAILURE);
             default: usage();
                 exit(EXIT_FAILURE);
        }
    }

    if(counter < 4) {
      G4cout << G4endl <<
      "**************************************************************************"
      << G4endl <<  
      "*                        Missing required parameters                     *"           
      << G4endl <<  
      "**************************************************************************"
      << G4endl << G4endl;  
      usage(); 
      exit(EXIT_FAILURE);
    }
}
