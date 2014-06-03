#ifndef MyEventAction_h
#define MyEventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class MyRunAction;
class MySteppingAction;


class MyEventAction: public G4UserEventAction {

public:

  MyEventAction();
  ~MyEventAction();

  virtual void BeginOfEventAction( const G4Event* evt );    
  virtual void EndOfEventAction( const G4Event* evt );    
  virtual void EndOfEventActionMaster( const G4Event* evt );    
  //01.25.2009 Xin Dong: The member field to store the result
  //of the event processed.
  G4double totalEdepAllParticles;

private:

  void instanciateRunAction();
  void instanciateSteppingAction();

  MyRunAction* theRunAction;
  MySteppingAction* theSteppingAction;

};

#endif
