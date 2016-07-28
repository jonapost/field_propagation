#ifndef PhysicsList_HH
#define PhysicsList_HH

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class PhysicsList: public G4VUserPhysicsList
{
public:
  PhysicsList();
  virtual ~PhysicsList();

  // Construct particles
  virtual void ConstructParticle();

  // Construct processes and register them
  virtual void ConstructProcess();

private:

};
#endif
