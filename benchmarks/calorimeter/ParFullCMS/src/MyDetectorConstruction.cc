#include "globals.hh"
#include "G4VisAttributes.hh"
#include "MyDetectorConstruction.hh"
#include "G4BooleanSolid.hh"
#include "G4CSGSolid.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "MyDetectorMessenger.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

// 04-June-2007 : Added to test Gabriele's tolerance.
#include "G4GeometryTolerance.hh"

//pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;

//01.25.2009 Xin Dong: Threads do not share this member data.
G4ThreadLocal G4FieldManager* MyDetectorConstruction::fieldMgr = 0;

//01.25.2009 Xin Dong: Threads do not share this member data.
G4ThreadLocal G4UniformMagField* MyDetectorConstruction::uniformMagField = 0;

//01.25.2009 Xin Dong: Threads do not share this member data.
G4ThreadLocal MyDetectorMessenger* MyDetectorConstruction::detectorMessenger = 0;

MyDetectorConstruction::MyDetectorConstruction() {

  ///testing GDML writer
  ///<<<<<<< MyDetectorConstruction.cc
  //parser.Read( "cms_out.gdml" );  //***LOOKHERE***
  ///parser.Write("cms_out.gdml",parser.GetWorldVolume()->GetLogicalVolume());
  //=======

  //  pthread_mutex_lock(&mut);

  parser = new G4GDMLParser;
  parser->Read( "cms.gdml" );   //***LOOKHERE***

  //  pthread_mutex_unlock(&mut);
  ///>>>>>>> 1.3
  detectorMessenger = new MyDetectorMessenger( this );
}

//01.25.2009 Xin Dong: Used by worker threads to achieve the partial
//effect similar to the member function Construct() invoked by the
//master thread.
void MyDetectorConstruction::SlaveMyDetectorConstruction()
{
  detectorMessenger = new MyDetectorMessenger(this);
}

//01.25.2009 Xin Dong: Use by worker threads to achieve the partial
//effect similar to the destructor invoked by the master thread.
void MyDetectorConstruction::SlaveDestroy()
{
  if (!uniformMagField)
    delete uniformMagField;
  if (detectorMessenger)
    delete detectorMessenger;
}

MyDetectorConstruction::~MyDetectorConstruction() {
  delete uniformMagField;
  delete detectorMessenger;
}

//01.25.2009 Xin Dong: Used by worker threads to achieve the partial
//effect similar to the member function Construct() invoked by the
//master thread.
G4VPhysicalVolume* MyDetectorConstruction::ConstructSlave() {
  fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  return fWorld;
}

void MyDetectorConstruction::ConstructSDandField() {
  //TODO: move from here 
  //TLS initialization if needed
  if ( !detectorMessenger ) detectorMessenger = new MyDetectorMessenger(this);
  if (! fieldMgr ) fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
}

G4VPhysicalVolume* MyDetectorConstruction::Construct() { 

  fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
 
  fWorld = (G4VPhysicalVolume *) parser->GetWorldVolume();

  fWorld->GetLogicalVolume()->SetVisAttributes (G4VisAttributes::Invisible);
  
  if( fWorld == 0 ) {
    G4ExceptionDescription msg;
    msg <<"World volume not set properly check your setup selection criteria or GDML input!";
    G4Exception("MyDetectorConstruction","ParFullCMS_001",FatalException,msg);
  }

  // ------------------------------------------------------------
  // 04-June-2007 : Added to test Gabriele's tolerance.

  // Calculate the dimensions of the World Volume by randomly
  // drawing points on its surface. 
   
  G4double minX = 0.0, maxX = 0.0;
  G4double minY = 0.0, maxY = 0.0;
  G4double minZ = 0.0, maxZ = 0.0;
  for ( int i = 0 ; i < 1000 ; i++ ) {
    G4ThreeVector surfacePoint = 
      fWorld->GetLogicalVolume()->GetSolid()->GetPointOnSurface();
    if ( surfacePoint.x() < minX ) minX = surfacePoint.x();
    if ( surfacePoint.x() > maxX ) maxX = surfacePoint.x();
    if ( surfacePoint.y() < minY ) minY = surfacePoint.y();
    if ( surfacePoint.y() > maxY ) maxY = surfacePoint.y();
    if ( surfacePoint.z() < minZ ) minZ = surfacePoint.z();
    if ( surfacePoint.z() > maxZ ) maxZ = surfacePoint.z();
  }

  // Consider the maximum extent of the World Volume.
  G4double maxWorldExtent = maxX - minX;
  if ( maxWorldExtent < ( maxY - minY ) ) maxWorldExtent = maxY - minY;
  if ( maxWorldExtent < ( maxZ - minZ ) ) maxWorldExtent = maxZ - minZ;

  G4cout << G4endl 
         << " WORLD VOLUME = " 
         << fWorld->GetName() << "\t"
         << fWorld->GetLogicalVolume()->GetSolid()->GetName() << G4endl
	 << " \t x-axis : " << minX/m << " , " << maxX/m << " metres" << G4endl
	 << " \t y-axis : " << minY/m << " , " << maxY/m << " metres" << G4endl
	 << " \t z-axis : " << minZ/m << " , " << maxZ/m << " metres" << G4endl
         << " \t -->  maxWorldExtent = " << maxWorldExtent/m << " metres" << G4endl 
         << G4endl;

  G4cout << G4endl
	 << " Compute tolerance = "
	 << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/m 
	 << " metres" << G4endl
	 << G4endl;

  // ------------------------------------------------------------

  delete parser;

  return fWorld;
}


void MyDetectorConstruction::SetMagField( const G4double fieldValue ) {
  if ( uniformMagField ) {
    delete uniformMagField;
  }
  if ( std::abs( fieldValue ) > 0.0 ) {
    // Apply a global uniform magnetic field along the Z axis.
    // Notice that only if the magnetic field is not zero, the Geant4
    // transportion in field gets activated.

    uniformMagField = new G4UniformMagField( G4ThreeVector( 0.0, 0.0, fieldValue ) );

    fieldMgr->SetDetectorField( uniformMagField );
    fieldMgr->CreateChordFinder( uniformMagField );

    G4cout << G4endl
           << " *** SETTING MAGNETIC FIELD : fieldValue = " << fieldValue / tesla
           << " Tesla *** " << G4endl 
	   << G4endl;

  } 
}
