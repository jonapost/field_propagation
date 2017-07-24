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

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
//#include "CLHEP/Units/SystemOfUnits.h"

#include "G4Step.hh"
#include "G4RunManager.hh"

#include "RunAction.hh"

#include "G4Track.hh"

double SteppingAction::fgTrackingCutInEnergy = 0.001; // tracking cut in GeV default value

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
    : fDetector((DetectorConstruction *)G4RunManager::GetRunManager()->GetUserDetectorConstruction()),
      fEventaction((EventAction *)G4RunManager::GetRunManager()->GetUserEventAction()) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *aStep) {

  if (!RunAction::fgScoreTypeFlag)
    return;

  unsigned long numPysLimSteps = 0;
  unsigned long numSecs = 0;
  unsigned long numTotal = 0;

  if (aStep->GetPreStepPoint()->GetKineticEnergy() > fgTrackingCutInEnergy) // get good number of secondaries
  {
    // determine which process happend
    G4VProcess const *g4proc = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    G4int g4procCode;
    G4String pNameStr;
    if (g4proc) // if NULL -> it was user special cut along the step that was limited by NOT transportation
    {
      g4procCode = g4proc->GetProcessType() * 1000 + g4proc->GetProcessSubType();
      pNameStr = g4proc->GetProcessName();
    } else {
      g4procCode = 0;
      pNameStr = "Nothing"; // the important is that the name is Not UserSpecialCut
    }
    if (RunAction::isTabPhys) { // running with G4 phys. list. convert G4ProcName to GV
      if (g4procCode != 1091) {
        ++numPysLimSteps;
        for (unsigned int i = 0; i < aStep->GetSecondaryInCurrentStep()->size(); ++i)
          if ((*aStep->GetSecondaryInCurrentStep())[i]->GetKineticEnergy() > fgTrackingCutInEnergy)
            ++numSecs;
      }
      if (RunAction::fgScoreTypeFlag > 1)
        fEventaction->FillHistSteps(aStep);
    } else {                                                              // G4 Physics
      if (g4procCode != 1091 && G4String("UserSpecialCut") != pNameStr) { // not transportation nt trackingCut
        ++numPysLimSteps;
        for (unsigned int i = 0; i < aStep->GetSecondaryInCurrentStep()->size(); ++i)
          if ((*aStep->GetSecondaryInCurrentStep())[i]->GetKineticEnergy() > fgTrackingCutInEnergy)
            ++numSecs;
      }
      if (RunAction::fgScoreTypeFlag > 1)
        fEventaction->FillHistSteps(aStep);
    }
    numTotal = 1;
  }
  // printf("FillPerStep: numPhys=%ld numSecs=%ld, numTotal=%ld\n", numPysLimSteps, numSecs, numTotal);
  if (RunAction::fgScoreTypeFlag > 0)
    fEventaction->FillPerSteps(numPysLimSteps, numSecs, numTotal, 1);

  // example of saving random number seed of this event, under condition
  //// if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();
}
