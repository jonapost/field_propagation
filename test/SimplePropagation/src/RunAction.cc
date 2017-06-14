#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4ParameterManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"

RunAction::RunAction(): G4UserRunAction()
{  
    G4RunManager* runManager = G4RunManager::GetRunManager();
    const G4VUserDetectorConstruction* detector = runManager->GetUserDetectorConstruction();

    pField = static_cast<const DetectorConstruction*>(detector)->GetField();
}



RunAction::~RunAction()
{}



void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
}



void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
    G4cout<<"Number of calls to field: "<<pField->GetCountCalls()<<G4endl;
}




