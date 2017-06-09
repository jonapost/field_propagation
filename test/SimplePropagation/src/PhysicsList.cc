#include "PhysicsList.hh"
#include "G4Proton.hh"

PhysicsList::PhysicsList(): G4VUserPhysicsList()
{
}


PhysicsList::~PhysicsList()
{
}


void PhysicsList::ConstructParticle()
{
  G4Proton::ProtonDefinition();
}


void PhysicsList::ConstructProcess()
{
  AddTransportation();
}
