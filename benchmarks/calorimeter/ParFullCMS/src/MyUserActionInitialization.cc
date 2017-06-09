#include "MyUserActionInitialization.hh"
#include "MyPrimaryGeneratorAction.hh"
#include "MyEventAction.hh"

void MyUserActionInitialization::Build() const {
  SetUserAction( new MyPrimaryGeneratorAction ); 
  SetUserAction( new MyEventAction ); 
}

void MyUserActionInitialization::BuildForMaster() const {
  //TODO: put here global run-manager with reductions, to be done
}
