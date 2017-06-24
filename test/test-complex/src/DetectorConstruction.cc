
#include "DetectorConstruction.hh"
#include "TGeoManager.h"
#include "TabulatedDataManager.hh"
#include "MaterialConverter.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"

#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4CashKarpRKF45.hh"

#include "G4TransportationManager.hh"

#include "CMSMagneticFieldG4.hh"
#include "G4CachedMagneticField.hh"

TGeoManager * DetectorConstruction::fgGeomMgrRoot= 0 ; // Pointer to the geometry manager   
CMSMagneticFieldG4* DetectorConstruction::fpMasterCmsField= 0;

DetectorConstruction::DetectorConstruction(G4VPhysicalVolume *setWorld, bool uniformMagField ) :
   fWorld(setWorld), fUseUniformField( uniformMagField ), fpMagField(0),
   fFieldFilename("cmsmagfield2015.txt"),
   fAllocatedStepper(nullptr),
   fMinStepField( 1.0e-2 * CLHEP::millimeter ),
   fDistanceConst( 1.0 * CLHEP::millimeter )
{   
   // fWorld = setWorld;
/*      char* gdmlFileName = getenv("VP_GEOM_GDML");
      if( gdmlFileName ){
        std::cout << " Creating empty TGeoManager by reading Root geometry from file " << gdmlFileName  << G4endl;
        fgGeomMgrRoot = TGeoManager::Import(gdmlFileName);
      } else {
        std::cout << " Creating empty TGeoManager " << std::endl;
        fgGeomMgrRoot = new TGeoManager();
      }
*/
//      std::cout << " Creating empty TGeoManager " << std::endl;
  fgGeomMgrRoot = new TGeoManager();

  TabulatedDataManager::SetTGeomManager( fgGeomMgrRoot );
  MaterialConverter::SetTGeomManager( fgGeomMgrRoot );
  MaterialConverter::Instance();// call just to initialize

  if( fUseUniformField ) { 
     // initialize magnetic field :: same value as in the prototype
     SetUniformBzMagField( fBzFieldValue );
  } else { 
     // To use the CMS field map - must create it from the field-dump file
     if( ! fpMasterCmsField ) 
        fpMasterCmsField= new CMSMagneticFieldG4( fFieldFilename );
     fpMagField= fpMasterCmsField;     
  }
}

DetectorConstruction::~DetectorConstruction()
{
  if( fpMagField != fpMasterCmsField ) delete fpMagField;
  delete fpMasterCmsField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetUniformBzMagField(G4double BzValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

//  if(fpMagField) delete fpMagField;                //delete the existing magn field

  fBzFieldValue=    BzValue;
  fUseUniformField= true;
  
  delete fpMagField;
  delete fieldMgr->GetChordFinder();
  fieldMgr->SetChordFinder(nullptr);
     
  if(BzValue!=0.)                        // create a new one if non nul
  {
    fpMagField = new G4UniformMagField(G4ThreeVector(0.,0.,BzValue));
    fieldMgr->SetDetectorField(fpMagField);
    // fieldMgr->CreateChordFinder(fpMagField);
  } else {
    fpMagField = 0;
    fieldMgr->SetDetectorField(fpMagField);
  }

  /***   The new way of creating it -- following examples/basic 2015
  G4ThreeVector fieldVector = G4ThreeVector(0.,0.,BzValue));
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(BzValue);
  fMagFieldMessenger->SetVerboseLevel(1);  
  ***/
}

void DetectorConstruction::CreateAndRegisterCMSfield()
{
   if( fpMagField != fpMasterCmsField ) delete fpMagField;

   if( ! fpMasterCmsField ) {
      std::cout << " DetectorConstruction::CreateAndRegisterCMSfield : initialising CMSMagneticFieldG4" << std::endl;
      fpMasterCmsField= new CMSMagneticFieldG4( fFieldFilename );
   }

   // fpMagField= fpMasterCmsField; // ->CloneOrSafeSelf();
   // Re-use the same object in all threads -- methods must remain 'truly' const

   fpMagField= new G4CachedMagneticField( fpMasterCmsField, fDistanceConst );
   // Field will not be re-evaluated if position moves by less than fDistanceConst
   // Note: The cache object will be different in each thread, as needed
   
   // Check the values
   double position[3] = { 0. , 0., 0. };
   double BfieldVal[3];
   fpMagField->GetFieldValue( position, BfieldVal ); 
   G4cout << " Test value of B field at position (0, 0, 0):  B(x,y,z) = ( "
          << (BfieldVal[0] / CLHEP::tesla) << " , " << (BfieldVal[1] / CLHEP::tesla) << " , "
          << (BfieldVal[2] / CLHEP::tesla)
          << " )  Tesla " << G4endl;
   
   G4FieldManager* fieldMgr
      = G4TransportationManager::GetTransportationManager()->GetFieldManager();
   // delete fieldMgr->GetChordFinder();  --> Done in CreateChordFinder already

   fieldMgr->SetDetectorField(fpMagField);
   // fieldMgr->CreateChordFinder(fpMagField);
}

G4ChordFinder*
DetectorConstruction::
CreateChordFinder( G4MagneticField*        magField,
                   G4MagIntegratorStepper* stepper )
{
  G4MagInt_Driver*        intgrDriver;
   
  //  Construct the Chord Finder
  //  by creating in inverse order the  Driver, the Stepper and EqRhs ...

  G4Mag_EqRhs *pEquation = new G4Mag_UsualEqRhs(magField);
  //  fEquation = pEquation;                            
  // fLastStepEstimate_Unconstrained = DBL_MAX;          // Should move q, p to
                                                     //    G4FieldTrack ??

  //  SetFractions_Last_Next( fFractionLast, fFractionNextEstimate);  
    // check the values and set the other parameters

  if( stepper == 0 )
  { 
     // stepper = new G4ClassicalRK4(pEquation);
     // G4cout << "DetectorConstruction:  Created Stepper ClassicalRK4" << G4endl;
     stepper = new G4CashKarpRKF45(pEquation);
     G4cout << "DetectorConstruction:  Created Stepper CashKarpRKF45" << G4endl;
     fAllocatedStepper= stepper;
  }
  else
  {
     G4cout << "DetectorConstruction:  Using pre-existing Stepper" << G4endl;     
     fAllocatedStepper= nullptr;
  }

  intgrDriver = new G4MagInt_Driver( fMinStepField,
                                     stepper, 
                                     stepper->GetNumberOfVariables() );

  return new G4ChordFinder( intgrDriver );
}


// Called in each thread to create a magnetic field
//
void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

   
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4FieldManager* fieldMgr
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();  
   
  if( fUseUniformField ){
     SetUniformBzMagField( fBzFieldValue );
  } else {
     CreateAndRegisterCMSfield();
  }

  fieldMgr->SetDetectorField(fpMagField);    
  
  if( fpMagField ){
     // fieldMgr->CreateChordFinder(fpMagField);  //  ClassicalRK4
     
     G4ChordFinder* cf= CreateChordFinder(fpMagField);  // Use Cash Karp
     fieldMgr->SetChordFinder(cf);
  }
  
  // G4AutoDelete::Register(fMagFieldMessenger);
}
