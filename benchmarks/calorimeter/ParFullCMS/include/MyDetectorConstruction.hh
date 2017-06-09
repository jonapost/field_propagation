#ifndef MyDetectorConstruction_H
#define MyDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"

#include "G4GDMLParser.hh"

#include "G4Threading.hh"

class G4VPhysicalVolume;
class G4FieldManager;
class G4UniformMagField;
class G4Material;
class MyDetectorMessenger;


class MyDetectorConstruction : public G4VUserDetectorConstruction {

public:

  MyDetectorConstruction();
  ~MyDetectorConstruction();

  G4VPhysicalVolume* Construct();

  void ConstructSDandField();

  //01.25.2009 Xin Dong: Used by worker threads to achieve the partial
  //effect similar to the member function Construct() invoked by the
  //master thread.
  G4VPhysicalVolume* ConstructSlave();

  void SetMagField( const G4double fieldValue );

  //01.25.2009 Xin Dong: Used by worker threads to achieve the partial
  //effect similar to the constructor implicitly invoked by the master
  //thread.
  void SlaveMyDetectorConstruction();

  //01.25.2009 Xin Dong: Use by worker threads to achieve the partial
  //effect similar to the destructor invoked by the master thread.
  void SlaveDestroy();

  
private:
  
  G4VPhysicalVolume* fWorld;

  G4GDMLParser* parser;

  //01.25.2009 Xin Dong: Threads do not share this member data.
  static G4ThreadLocal G4FieldManager* fieldMgr;
  // Pointer to the field manager.

  //01.25.2009 Xin Dong: Threads do not share this member data.
  static G4ThreadLocal G4UniformMagField* uniformMagField; 
  // Pointer to the uniform magnetic field.
  
  //01.25.2009 Xin Dong: Threads do not share this member data.
  static G4ThreadLocal MyDetectorMessenger* detectorMessenger;
  // Pointer to the Messenger.

};

#endif

