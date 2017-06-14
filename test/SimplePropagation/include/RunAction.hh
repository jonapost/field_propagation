#ifndef RunAction_HH
#define RunAction_HH

#include "G4UserRunAction.hh"
#include "G4Parameter.hh"
#include "globals.hh"
#include "G4Run.hh"
#include "G4CachedMagneticField.hh"


class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void AddEdep (G4double edep); 

  private:

    G4CachedMagneticField* pField;
};

#endif

