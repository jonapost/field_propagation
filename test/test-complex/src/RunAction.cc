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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "TabulatedDataManager.hh"
#include "G4Step.hh"
#include "CMSApp.hh"

#ifdef USE_VECGEOM_NAVIGATOR
#undef USE_VECGEOM_NAVIGATOR
#define RESTORE_USE_VECGEOM_NAVIGATOR
#endif

#include "TPartIndex.h"

#ifdef RESTORE_USE_VECGEOM_NAVIGATOR
#define USE_VECGEOM_NAVIGATOR
#endif

#include <time.h>



G4bool RunAction::isTabPhys = FALSE;     // running with TABPHYS ?
G4int  RunAction::fgScoreTypeFlag = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction():
  fRunTime(0),
  fNumPysLimStepsRun(0),   // number of steps limited by physics
  fNumPrimsRun(0),         // number of primaries  
  fNumSecsRun(0),          // number of secondaries
  fNumTotalStepsRun(0),    // total number of steps   
  fNumAllStepsRun(0),      // number of ALL steps
  // fSumTime(0),
  fCMSApp(0)
{
  fRunTime=clock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  if(fgScoreTypeFlag > 1){ 
    fCMSApp = new CMSApp();
    fCMSApp->Initialize();
  }

  G4int sevent = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
  std::cerr << "<----- Run " << aRun->GetRunID() << " has started and will simulate " <<
            sevent << " events ! ----->\n";

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  fNumPysLimStepsRun  = 0;   // number of steps limited by physics
  fNumPrimsRun        = 0;   // number of primaries  
  fNumSecsRun         = 0;   // number of secondaries
  fNumTotalStepsRun   = 0;   // total number of steps   
  fNumAllStepsRun     = 0;   // number of ALL steps

  clock_t startClock= clock();
  G4cout << " Run started - " << ((G4double) (startClock - fRunTime)) / CLOCKS_PER_SEC
         << " [sec] since RunAction was created." << G4endl;
  
  fRunTime= clock(); // startClock;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillPerEvent(unsigned long nphyssteps, unsigned long nprims, 
                             unsigned long nsecs, unsigned long ntotal, unsigned long nall){

  fNumPysLimStepsRun += nphyssteps;  // number of steps limited by physics
  fNumPrimsRun       += nprims;      // number of primaries  
  fNumSecsRun        += nsecs;       // number of secondaries
  fNumTotalStepsRun  += ntotal;      // total number of steps   
  fNumAllStepsRun    += nall;        // number of ALL steps

}


void RunAction::FillHist(const G4Step *step) {fCMSApp->SteppingAction(step);}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  std::cerr<<"\n";
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
 
  fRunTime= clock() - fRunTime;
  std::cout<<"\n-------------------- Run time in [s] --------------------------\n"
	   <<"  Total run time= " << ((G4double)(clock()-fRunTime))/CLOCKS_PER_SEC 
            <<"\n--------------------------------------------------------------\n"
            << std::endl;

  if(fgScoreTypeFlag > 0) {
    std::cout<<"\n-----------------------------------------------------------------"        <<std::endl;
    std::cout<<"======  TOTAL NUMBER OF PRIMARY TRANSPORTED      = " << fNumPrimsRun        << std::endl;
    std::cout<<"======  TOTAL NUMBER OF SECONDARY TRANSPORTED    = " << fNumSecsRun         << std::endl;
    std::cout<<"======  TOTAL NUMBER OF STEPS                    = " << fNumTotalStepsRun   << std::endl;
    std::cout<<"======  TOTAL NUMBER OF STEPS LIMITED BY PHYSICS = " << fNumPysLimStepsRun  << std::endl;
    //std::cout<<"======  TOTAL NUMBER OF ALL STEPS                = " << fNumAllStepsRun     << std::endl;
    std::cout<<"-----------------------------------------------------------------"          <<std::endl;
  }

  if(fgScoreTypeFlag > 1)
    fCMSApp->EndOfRunAction(fNumPrimsRun);

  G4int sevent = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
  std::cerr << "<----- Run " << aRun->GetRunID() << " has finished the simulation of " << sevent << " events! ---->\n";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

