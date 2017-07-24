#ifndef EventAction_HH
#define EventAction_HH

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "RunAction.hh"

/// Event action class
///

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);


  private:
};


#endif

    
