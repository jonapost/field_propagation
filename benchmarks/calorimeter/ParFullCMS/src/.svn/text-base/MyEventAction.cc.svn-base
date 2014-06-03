#include "MyEventAction.hh"
#include "G4ios.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "MyRunAction.hh"
#include "MySteppingAction.hh"
#include "G4SystemOfUnits.hh"


MyEventAction::MyEventAction() : 
  theRunAction( 0 ), theSteppingAction( 0 ) {
  instanciateRunAction();
  instanciateSteppingAction();
}


MyEventAction::~MyEventAction() {}


void MyEventAction::BeginOfEventAction( const G4Event* ) { 
  //G4cout << "\n---> Begin of event: " << evt->GetEventID() << G4endl;
}


void MyEventAction::EndOfEventAction( const G4Event* evt ) {

  // Get from the stepping action the total deposited energy in
  // the event in the whole experimental hall.
  //01.25.2009 Xin Dong: Assign the member field with the result
  //of the event processed. This result will then be sent to the
  //master thread.
  totalEdepAllParticles = theSteppingAction->getTotalEdepAllParticles();

  G4cout << " ---  MyEventAction::EndOfEventAction  ---    event = " 
	 << evt->GetEventID() << "  Edep=" << totalEdepAllParticles/GeV
	 << " GeV" << G4endl;

  //  theRunAction->updateEndOfEvent( totalEdepAllParticles );
  theSteppingAction->reset();

}

void MyEventAction::EndOfEventActionMaster( const G4Event* ) {
  theRunAction->updateEndOfEvent( totalEdepAllParticles );
}


void MyEventAction::instanciateRunAction() {
  G4UserRunAction* theUserAction = const_cast< G4UserRunAction* >
    ( G4RunManager::GetRunManager()->GetUserRunAction() );
  if ( theUserAction == 0 ) {
    theRunAction = new MyRunAction;
    G4RunManager::GetRunManager()->SetUserAction( theRunAction );
  }
}


void MyEventAction::instanciateSteppingAction() {
  G4UserSteppingAction* theUserAction = const_cast< G4UserSteppingAction* >
    ( G4RunManager::GetRunManager()->GetUserSteppingAction() );
  if ( theUserAction == 0 ) {
    theSteppingAction = new MySteppingAction;
    G4RunManager::GetRunManager()->SetUserAction( theSteppingAction );
  }
}
