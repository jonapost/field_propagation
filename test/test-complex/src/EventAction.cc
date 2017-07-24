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

#include "EventAction.hh"

#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "TabulatedDataManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction():
   fRunAct((RunAction*)G4RunManager::GetRunManager()->GetUserRunAction())
{
  fNumPysLimStepsEvent  = 0;   // number of steps limited by physics
  fNumSecsEvent         = 0;   // number of secondaries
  fNumPrimsEvent        = 0;   // number of primaries
  fNumTotalStepsEvent   = 0;   // total number of steps (possible limit condition)
  fNumAllStepsEvent     = 0;   // number of all steps 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
 std::cout<< "----- Event number = "        << evt->GetEventID() 
          << "  with number of primaries = "<< evt->GetPrimaryVertex()->GetNumberOfParticle()
          << std::endl;

  fNumPysLimStepsEvent  = 0;   // number of steps limited by physics
  fNumSecsEvent         = 0;   // number of secondaries
  fNumTotalStepsEvent   = 0;   // number of steps
  fNumAllStepsEvent     = 0;

  fNumPrimsEvent        = evt->GetPrimaryVertex()->GetNumberOfParticle();   // number of secondaries  
}

void EventAction::FillPerSteps(unsigned long nphyssteps, unsigned long nsecs, unsigned long ntotal,
                               unsigned long nall){
  fNumPysLimStepsEvent += nphyssteps;   // number of steps limited by physics
  fNumSecsEvent        += nsecs;        // number of secondaries
  fNumTotalStepsEvent  += ntotal;
  fNumAllStepsEvent    += nall;
}

void EventAction::EndOfEventAction(const G4Event* /*evt*/){
  fRunAct->FillPerEvent(fNumPysLimStepsEvent, fNumPrimsEvent, fNumSecsEvent, fNumTotalStepsEvent, fNumAllStepsEvent);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::FillHistSteps(const G4Step *step){ fRunAct->FillHist(step); }

