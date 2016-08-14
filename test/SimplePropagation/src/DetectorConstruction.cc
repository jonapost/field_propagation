#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"


//field
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
//#include "G4BogackiShampine45.hh"
#include "G4PropagatorInField.hh"
//#include "G4DormandPrinceRK78.hh"
#include "G4SimpleLocator.hh"
#include "G4BSChordFinder.hh"
#include "G4RKChordFinder.hh"
#include "G4MagIntegratorDriver.hh"
#include "BulirschStoerDriver.hh"

#include "BulirschStoerDenseDriver.hh"
#include "G4BogackiShampine45DenseDriver.hh"
#include "G4RevisedChordFinder.hh"

using namespace CLHEP;

DetectorConstruction::DetectorConstruction():
    G4VUserDetectorConstruction(),
    fMinChordStep(0.0001)
{
    fpChordFinder = nullptr;
    field = new G4UniformMagField(G4ThreeVector(0.*tesla, 1.*tesla, 0.*tesla));
    fpField = new G4CachedMagneticField(field, 0);
}



DetectorConstruction::~DetectorConstruction()
{
    delete fpChordFinder;
    delete field;
    delete fpField;
}


void DetectorConstruction::ConstructField()
{
    G4Mag_UsualEqRhs *pEquation = new G4Mag_UsualEqRhs(fpField);
    G4FieldManager *globalFieldManager =
            G4TransportationManager::GetTransportationManager()->GetFieldManager();


    G4PropagatorInField * globalPropagatorInField =
            G4TransportationManager::GetTransportationManager()->GetPropagatorInField();

    G4Navigator* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    globalPropagatorInField->SetIntersectionLocator(new G4SimpleLocator(navigator));

    G4cout<<globalFieldManager->GetDeltaOneStep()<<"   "<<globalFieldManager->GetMinimumEpsilonStep()<<"   "<<globalFieldManager->GetDeltaIntersection()<<G4endl;

    globalFieldManager->SetDeltaOneStep(globalFieldManager->GetDeltaOneStep()/1e2);
    globalFieldManager->SetMinimumEpsilonStep(globalFieldManager->GetMinimumEpsilonStep()/1e2);
    globalFieldManager->SetDeltaIntersection(globalFieldManager->GetDeltaIntersection()/1e2);

    G4cout<<"DeltaOneStep: "<<globalFieldManager->GetDeltaOneStep()
          <<"  DeltaIntersection: "<<globalFieldManager->GetDeltaIntersection()<<G4endl;

    globalPropagatorInField->SetMaxLoopCount( 10000 );
    G4cout
      << "PropagatorInField parameter(s) are: " << G4endl
      << " SetMaxLoopCount=" << globalPropagatorInField->GetMaxLoopCount()
      << " minEpsilonStep= " << globalPropagatorInField->GetMinimumEpsilonStep() << " "
      << " maxEpsilonStep= " << globalPropagatorInField->GetMaximumEpsilonStep() << " "
      << G4endl;

    globalFieldManager->SetDetectorField(fpField);

    // globalFieldManager->SetMinimumEpsilonStep( 5.0e-7 );    // Old value
    // globalFieldManager->SetMaximumEpsilonStep( 0.05 );      // FIX - old value
    // globalFieldManager->SetDeltaOneStep( 0.25 * mm );       // original value
    // globalFieldManager->SetDeltaIntersection( 0.10 * mm );  // original value

    G4cout << "Field Manager's parameters are "
       << " minEpsilonStep= " << globalFieldManager->GetMinimumEpsilonStep() << " "
       << " maxEpsilonStep= " << globalFieldManager->GetMaximumEpsilonStep() << " "
       << " deltaOneStep=   " << globalFieldManager->GetDeltaOneStep() << " "
       << " deltaIntersection= " << globalFieldManager->GetDeltaIntersection()
       << G4endl;

    //G4MagIntegratorStepper* pStepper = new G4ClassicalRK4(pEquation);
    //G4MagIntegratorStepper* pStepper = new G4CashKarpRKF45(pEquation);
    //G4MagIntegratorStepper* pStepper = new G4BogackiShampine45(pEquation);
    //G4MagIntegratorStepper* pStepper = new G4DormandPrinceRK78(pEquation);
    //G4VIntegrationDriver* pDriver = new G4MagInt_Driver(fMinChordStep, pStepper);
    G4VIntegrationDriver* pDriver = new BulirschStoerDenseDriver(fMinChordStep, pEquation);
    fpChordFinder = new G4RevisedChordFinder(pDriver);
    fpChordFinder->SetVerbose(1);
    G4cout<<"DeltaChord: "<<fpChordFinder->GetDeltaChord()<<G4endl;
    fpChordFinder->SetDeltaChord(fpChordFinder->GetDeltaChord());
    globalFieldManager->SetChordFinder(fpChordFinder);

}


G4VPhysicalVolume* DetectorConstruction::Construct()
{  
    ConstructField();

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;


  // World
  G4double world_size = 200*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                                      //its name
       0.5*world_size, 0.5*world_size, 0.5*world_size);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  // Cube
  G4double cube_size = 20*cm;
  G4Material* cube_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidCube =
    new G4Box("Cube",                    //its name
        cube_size, 0.5*cube_size, cube_size); //its size
      
  G4LogicalVolume* logicCube =
    new G4LogicalVolume(solidCube,            //its solid
                        cube_mat,             //its material
                        "Cube");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicCube,                //its logical volume
                    "Cube",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 



  // Cube2
  //G4double cube_size = 20*cm;
  //G4Material* cube_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidCube2 =
    new G4Box("Cube2",                    //its name
        cube_size, 0.5*cube_size, cube_size); //its size

  G4LogicalVolume* logicCube2 =
    new G4LogicalVolume(solidCube2,            //its solid
                        cube_mat,             //its material
                        "Cube2");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, cube_size*2, 0),         //at (0,0,0)
                    logicCube2,                //its logical volume
                    "Cube2",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

                

  //always return the physical World
  return physWorld;
}


