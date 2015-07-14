//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: testPropagateMagField.cc 69604 2013-05-08 17:39:18Z japost $
//
//  
//
// Started from testG4Navigator1.cc,v 1.7 1996/08/29 15:42 pkent 
//   Locate & Step within simple boxlike geometry, both
//   with and without voxels. Parameterised volumes are included.

#include <assert.h>
#include <ErrorComputer.hh>
#include <stdio.h>
#include <time.h>
// #include "ApproxEqual.hh"

// Global defs
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4ios.hh"
#include <iomanip>

#include "Mag_UsualEqRhs_IntegrateByTime.hh"
#include "MagIntegratorStepperbyTime.hh"

#include <iostream>
#include <vector>

using namespace std;

#define BUFFER_MAX_LEN 100000

// Sample Parameterisation
class G4LinScale : public G4VPVParameterisation
{
	virtual void ComputeTransformation(const G4int n,
			G4VPhysicalVolume* pRep) const
	{
		pRep->SetTranslation(G4ThreeVector(0,(n-1)*15,0));
	}

	virtual void ComputeDimensions(G4Box &pBox,
			const G4int n,
			const G4VPhysicalVolume* ) const
	{
		pBox.SetXHalfLength(10);
		pBox.SetYHalfLength(5+n);
		pBox.SetZHalfLength(5+n);
	}

	virtual void ComputeDimensions(G4Tubs &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Trd &, 
			const G4int,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Cons &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Trap &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Hype &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Orb &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Sphere &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Torus &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Para &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Polycone &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
	virtual void ComputeDimensions(G4Polyhedra &,
			const G4int ,
			const G4VPhysicalVolume*) const {}
};
G4LinScale myParam;

// Build simple geometry:
// 4 small cubes + 1 slab (all G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{

	G4Box *myHugeBox=  new G4Box("huge box",15*m,15*m,25*m);
	G4Box *myBigBox=   new G4Box("big cube",10*m,10*m,10*m);
	G4Box *mySmallBox= new G4Box("smaller cube",2.5*m,2.5*m,2.5*m);
	G4Box *myTinyBox=  new G4Box("tiny  cube",.25*m,.25*m,.25*m);

	// G4Box *myVariableBox=
	new G4Box("Variable Box",10,5,5);

	//  World Volume
	//
	G4LogicalVolume *worldLog=new G4LogicalVolume(myHugeBox,0,
			"World",0,0,0);
	// Logical with no material,field,
	// sensitive detector or user limits

	G4PVPlacement *worldPhys=new 
		G4PVPlacement(0,G4ThreeVector(0,0,0), "World",worldLog,
				0,false,0);
	// Note: no mother pointer set

	//  Create the logical Volumes
	//
	//  G4LogicalVolume(*pSolid, *pMaterial, Name, *pField, *pSDetector, *pULimits)
	//
	G4LogicalVolume *BigBoxLog=new G4LogicalVolume(myBigBox,0,
			"Crystal Box (large)",0,0,0);
	G4LogicalVolume *smallBoxLog=new G4LogicalVolume(mySmallBox,0,
			"Crystal Box (small)");
	G4LogicalVolume *tinyBoxLog=new G4LogicalVolume(myTinyBox,0,
			"Crystal Box (tiny)");


	//  Place them.
	//
	//  1) Two big boxes in the world volume
	//
	// G4PVPlacement *BigTg1Phys=
	new G4PVPlacement(0,G4ThreeVector(0,0,-15*m),
			"Big Target 1",BigBoxLog,
			worldPhys,false,0);
	// G4PVPlacement *BigTg2Phys=
	new G4PVPlacement(0,G4ThreeVector(0,0, 15*m),
			"Big Target 2",BigBoxLog,
			worldPhys,false,0);

	//  2) Four (medium) boxes in X & Y near the origin of the world volume
	//
	// G4PVPlacement *MedTg3a_Phys=
	new G4PVPlacement(0,G4ThreeVector(0, 7.5*m,0),
			"Target 3a",smallBoxLog,
			worldPhys,false,0);
	// G4PVPlacement *MedTg3b_Phys=
	new G4PVPlacement(0,G4ThreeVector(0,-7.5*m,0),
			"Target 3b",smallBoxLog,
			worldPhys,false,0);
	// G4PVPlacement *MedTg3c_Phys=
	new G4PVPlacement(0,G4ThreeVector(-7.5*m,0,0),
			"Target 3c",smallBoxLog,
			worldPhys,false,0);
	// G4PVPlacement *MedTg3d_Phys=
	new G4PVPlacement(0,G4ThreeVector( 7.5*m,0,0),
			"Target 3d",smallBoxLog,
			worldPhys,false,0);


	//  3) Eight small boxes around the origin of the world volume 
	//        (in +-X, +-Y & +-Z)
	//
	// G4PVPlacement *SmTg4a_Phys=
	new G4PVPlacement
		(0,G4ThreeVector( 0.3*m, 0.3*m,0.3*m), "Target 4a",tinyBoxLog,
		 worldPhys,false,0);
	// G4PVPlacement *SmTg4b_Phys=
	new G4PVPlacement
		(0,G4ThreeVector( 0.3*m,-0.3*m,0.3*m), "Target 4b",tinyBoxLog,
		 worldPhys,false,0);
	// G4PVPlacement *SmTg4c_Phys=
	new G4PVPlacement
		(0,G4ThreeVector(-0.3*m,-0.3*m,0.3*m), "Target 4c",tinyBoxLog,
		 worldPhys,false,0);
	// G4PVPlacement *SmTg4d_Phys=
	new G4PVPlacement
		(0,G4ThreeVector(-0.3*m, 0.3*m,0.3*m), "Target 4d",tinyBoxLog,
		 worldPhys,false,0);

	// G4PVPlacement *SmTg4e_Phys=
	new G4PVPlacement
		(0,G4ThreeVector( 0.3*m, 0.3*m,-0.3*m), "Target 4e",tinyBoxLog,
		 worldPhys,false,0);
	// G4PVPlacement *SmTg4f_Phys=
	new G4PVPlacement
		(0,G4ThreeVector( 0.3*m,-0.3*m,-0.3*m), "Target 4f",tinyBoxLog,
		 worldPhys,false,0);
	// G4PVPlacement *SmTg4g_Phys=
	new G4PVPlacement
		(0,G4ThreeVector(-0.3*m,-0.3*m,-0.3*m), "Target 4g",tinyBoxLog,
		 worldPhys,false,0);
	// G4PVPlacement *SmTg4h_Phys=
	new G4PVPlacement
		(0,G4ThreeVector(-0.3*m, 0.3*m,-0.3*m), "Target 4h",tinyBoxLog,
		 worldPhys,false,0);

	return worldPhys;
}


#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4CachedMagneticField.hh"

#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ConstRK4.hh"
#include "G4NystromRK4.hh"
#include "G4HelixMixedStepper.hh"
#include "globals.hh"

#include "G4VCurvedTrajectoryFilter.hh"

#include "ChawlaSharmaRKNstepper.hh"
#include "FineRKNG34.hh"
#include "DormandPrince745.hh"

//=============test template mode================
/*
#include "TMagFieldEquation.hh"
#include "TCashKarpRKF45.hh"
#include "TCachedMagneticField.hh"
#include "TQuadrupoleMagField.hh"
#include "TClassicalRK4.hh"
#include "TSimpleHeum.hh"
#include "TSimpleRunge.hh"
#include "TExplicitEuler.hh"


//typedef G4CachedMagneticField Field_t;
//typedef TCachedMagneticField<G4QuadrupoleMagField> Field_t;
typedef TCachedMagneticField<TQuadrupoleMagField> Field_t;
typedef TMagFieldEquation<Field_t> Equation_t;
typedef TCashKarpRKF45<Equation_t, 6> Stepper_t;
typedef TClassicalRK4<Equation_t, 6> StepperRK4_t;
typedef TSimpleHeum<Equation_t, 6> StepperHeum_t;
typedef TSimpleRunge<Equation_t, 6> StepperRunge_t;
typedef TExplicitEuler<Equation_t, 6> StepperExEuler_t;


TQuadrupoleMagField   tQuadrupoleMagField( 10.*tesla/(50.*cm) ); 
//G4QuadrupoleMagField   tQuadrupoleMagField( 10.*tesla/(50.*cm) ); 
Field_t  tMagField( &tQuadrupoleMagField, 1.0 * cm); 
//===============================================
*/


//G4UniformMagField      uniformMagField( G4ThreeVector(0., 0., -1.*tesla) );
//G4CachedMagneticField  myMagField( &uniformMagField, 1.0 * cm);
//G4String   fieldName("Uniform -1.0 Tesla");

G4QuadrupoleMagField   quadrupoleMagField( 10.*tesla/(50.*cm) );
G4CachedMagneticField  myMagField( &quadrupoleMagField, 1.0 * cm);
G4String   fieldName("Cached Quadropole field, 20T/meter, cache=1cm");

G4FieldManager* SetupField(G4int type)
{
	G4FieldManager   *pFieldMgr;
	G4ChordFinder    *pChordFinder;

	Mag_UsualEqRhs_IntegrateByTime *fEquation = new Mag_UsualEqRhs_IntegrateByTime(&myMagField);
	G4MagIntegratorStepper *pStepper;

	//G4cout << " Setting up field of type: " << fieldName << G4endl;
	switch ( type ) 
	{
	   case -2: pStepper = new MagIntegratorStepper_byTime<DormandPrince745>(fEquation);
	      break;
	   case -1: pStepper = new MagIntegratorStepper_byTime<G4SimpleRunge>(fEquation);
         break;
	   case 0: pStepper = new MagIntegratorStepper_byTime<G4ClassicalRK4>(fEquation);
         break;
	   case 1: pStepper = new MagIntegratorStepper_byTime<G4CashKarpRKF45>(fEquation);
         break;
      case 2: pStepper = new MagIntegratorStepper_byTime<ChawlaSharmaRKNstepper>(fEquation);
         break;
      case 3: pStepper = new MagIntegratorStepper_byTime<FineRKNG34>(fEquation);
         break;
	   }

    pFieldMgr= G4TransportationManager::GetTransportationManager()->
       GetFieldManager();

    pFieldMgr->SetDetectorField( & myMagField );

    pChordFinder = new G4ChordFinder( & myMagField,
				      1.0e-2 * mm,
				      pStepper);


    pChordFinder->SetVerbose(0);  // ity();

    pFieldMgr->SetChordFinder( pChordFinder );

    return    pFieldMgr;
}

#include "G4SimpleLocator.hh"
#include "G4BrentLocator.hh"
#include "G4MultiLevelLocator.hh"

G4PropagatorInField*  SetupPropagator( G4int type)
{
    // G4FieldManager* fieldMgr=
    SetupField( type ) ;

    // G4ChordFinder  theChordFinder( &MagField, 0.05*mm ); // Default stepper
 
    G4PropagatorInField *thePropagator = 
      G4TransportationManager::GetTransportationManager()->
       GetPropagatorInField ();

    // Let us test the new Minimum Epsilon Step functionality
    // thePropagator -> SetMinimumEpsilonStep( 1.0e-3 ) ; 
    // thePropagator -> SetMaximumEpsilonStep( 1.0e-5 ) ; 

    G4Navigator *theNavigator= G4TransportationManager::GetTransportationManager()->
       GetNavigatorForTracking();
    // Test the options for Locator
    G4VIntersectionLocator *pLocator=0;
    //G4cout << "Over-riding  PropagatorInField to use ";
    pLocator= new G4MultiLevelLocator(theNavigator); //G4cout << "Multi"; // default
    // pLocator= new G4SimpleLocator(theNavigator); //G4cout << "Simple";
    // pLocator= new G4BrentLocator(theNavigator); //G4cout << " Brent "; 
    //G4cout << " Locator. ( In the unit test code. ) " << G4endl;

    thePropagator->SetIntersectionLocator(pLocator);

    return thePropagator;
}

G4PropagatorInField *pMagFieldPropagator=0; 

//
// Test Stepping
//
void testG4PropagatorInField(G4VPhysicalVolume*,     // *pTopNode,
			       G4int             	type,
				   G4double			 	step_distance,
				   G4double				step_no,
				   vector< vector<G4double> > &pos_mom_time_buffer ) {

    G4Navigator   *pNavig= G4TransportationManager::
                    GetTransportationManager()-> GetNavigatorForTracking();
    
    pMagFieldPropagator = SetupPropagator(type);

    G4double particleCharge= +1.0;  // in e+ units
    G4double spin=0.0;              // ignore the spin
    G4double magneticMoment= 0.0;   // ignore the magnetic moment

        // pNavig->SetWorldVolume(pTopNode);

    G4VPhysicalVolume *located;
    G4double step_len, physStep, safety;
    G4ThreeVector xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    
    // physStep=kInfinity;
    G4ThreeVector Position(0.,0.,0.); // changed here
    G4ThreeVector UnitMomentum(0.,0.6,0.8);  
    G4ThreeVector EndPosition, EndUnitMomentum;

//
// Test location & Step computation
//  
    /* assert(located->GetName()=="World"); */
    if( std::fabs(UnitMomentum.mag() - 1.0) > 1.e-8 ) 
    {
      G4cerr << "UnitMomentum.mag() - 1.0 = " << UnitMomentum.mag() - 1.0 <<
	G4endl;
    }

       physStep=  step_distance * mm ;  // millimeters

       Position = G4ThreeVector(0.,0.,0.);
       UnitMomentum = (G4ThreeVector(0.3,0.6,0.8)).unit();

       G4double momentum = (0.5) * proton_mass_c2;

       G4ThreeVector Momentum = UnitMomentum;
       Momentum *= momentum;

       G4double kineticEnergy =  momentum*momentum /
                  ( std::sqrt( momentum*momentum + proton_mass_c2 * proton_mass_c2 ) 
		    + proton_mass_c2 );
       G4double velocity = momentum / ( proton_mass_c2 + kineticEnergy );
       G4double labTof= 10.0*ns, properTof= 0.1*ns;
       G4ThreeVector Spin(1.0, 0.0, 0.0);

       G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                                     spin=0.0,
                                     magneticMoment=0.0);

         G4EquationOfMotion* equationOfMotion =
            ( pMagFieldPropagator->GetChordFinder()->GetIntegrationDriver()->GetStepper() )
            ->GetEquationOfMotion();
       // Momentum in Mev/c ?
       // pMagFieldPropagator
       equationOfMotion->SetChargeMomentumMass(
		      chargeState,//+1,                    // charge in e+ units
		      momentum, 
		      proton_mass_c2 ); // Changed here

       //G4double radius_curvature = ((proton_mass_c2 + kineticEnergy) / momentum)
       //                        * sqrt( y_initial[3]*y_initial[3] + y_initial[4]*y_initial[4] )
       //                        / ( particleCharge * 1. ); // In our case only y component is perpendicular to B field
                           // We don't multiply by mass because it is included in momentum

       //cout << "Radius of Curvature: " << radius_curvature << endl;

       //pMagFieldPropagator -> GetChordFinder() -> SetMass();
       //pMagFieldPropagator -> GetChordFinder() -> setup_output_buffer( Position,
       //                              Momentum, pos_mom_time_buffer);

       clock_t total = 0;
	   for( int istep=0; istep < step_no; istep++ ){
	  located= pNavig->LocateGlobalPointAndSetup(Position);

          G4FieldTrack  initTrack( Position, 
				   UnitMomentum,
				   0.0,            // starting S curve len
				   kineticEnergy,
				   proton_mass_c2,
				   velocity,
				   labTof, 
				   properTof,
				   0              // or &Spin
				   ); 
		  clock_t t;
		  t = clock(); 

	  step_len=pMagFieldPropagator->ComputeStep( initTrack, 
						     physStep, 
						     safety,
						     located);


	  //bufferA = pMagFieldPropagator -> GetChordFinder() -> GetBuffer();
	  //pMagFieldPropagator -> GetChordFinder() -> Reset_Buffer();


	  total += clock() - t;
	  //       --------------------
	  EndPosition=     pMagFieldPropagator->EndPosition();
	  EndUnitMomentum= pMagFieldPropagator->EndMomentumDir();
	  //       --------
	  
	  if( std::fabs(EndUnitMomentum.mag2() - 1.0) > 1.e-8 )
	    G4cerr << "EndUnitMomentum.mag2() - 1.0 = " <<
	      EndUnitMomentum.mag2() - 1.0 << G4endl;

	  G4ThreeVector MoveVec = EndPosition - Position;
	  assert( MoveVec.mag() < physStep*(1.+1.e-9) );

	  //report_endPV(Position, UnitMomentum, step_len, physStep, safety,
	  //             EndPosition, EndUnitMomentum, istep, located );

	  assert(safety>=0);
	  pNavig->SetGeometricallyLimitedStep();
	  // pMagFieldPropagator->SetGeometricallyLimitedStep();

	  Position= EndPosition;
	  UnitMomentum= EndUnitMomentum;
	  //if ( (istep + 1) % factor == 0)
	     //physStep *= 2;
	  //physStep *= 1.;
     } // ...........................  end for ( istep )
     //G4cout << "=============="<<total<<"================="<<G4endl;
	  //myMagField.ReportStatistics();

    //return pMagFieldPropagator -> GetChordFinder() -> GetCounter();
    //return(1);
}

// Main program
// -------------------------------
int main(int argc, char **argv)
{

    G4VPhysicalVolume *myTopNode;
    G4int type, step_distance, step_no, optim, optimSaf;
    G4bool optimiseVoxels=true;
    G4bool optimisePiFwithSafety=true;

    G4double step_distance_input;

    //G4cout << " Arguments:  stepper-no  optimise-Voxels optimise-PiF-with-safety" << G4endl;

    if( argc >= 2 ){
       type = atoi(argv[1]);
    }else{
		type = 3;
	}
    if( argc >= 3 ){
    	step_distance_input = atof(argv[2]);
    }else{
    	step_distance_input = 2.0;
    }
    if( argc >= 4 ){
        	step_no = atoi(argv[3]);
      }else{
      step_no = 3;
}


    vector< vector<G4double> > pos_mom_time_buffer;

   G4double err_pos_mag2, err_pos_mag;

   G4double momentum = 0.5 * proton_mass_c2;

   G4double kineticEnergy =  momentum*momentum /
          ( std::sqrt( momentum*momentum + proton_mass_c2 * proton_mass_c2 )
         + proton_mass_c2 );

   G4double radius_curvature;

   myTopNode=BuildGeometry();	// Build the geometry
 
    G4Navigator *pNavig= G4TransportationManager::
                    GetTransportationManager()-> GetNavigatorForTracking();
    pNavig->SetWorldVolume(myTopNode);

    G4GeometryManager::GetInstance()->CloseGeometry(false);

    // Setup the propagator (will be overwritten by testG4Propagator ...)
    pMagFieldPropagator = SetupPropagator(type);
    pMagFieldPropagator->SetUseSafetyForOptimization(optimisePiFwithSafety);

    testG4PropagatorInField( myTopNode, type,
          step_distance_input, step_no, pos_mom_time_buffer );

    G4GeometryManager::GetInstance()->OpenGeometry();


/*
    for (int i = 0; i < pos_mom_time_buffer.size(); i ++) {
       for (int j = 0; j < 3; j ++) {
          cout << pos_mom_time_buffer[i][j] << ",";
       }
       cout << pos_mom_time_buffer[i][9] << endl;
    }
*/
    /*
    for (int i = 0; i < counter2; i ++) {
       radius_curvature = (proton_mass_c2)
                                  * sqrt( buffer2[i][3]*buffer2[i][3] + buffer2[i][4]*buffer2[i][4] )
                                  / ( 1. * 1. );
       cout << radius_curvature << endl;
       for (int j = 0; j < 3; j ++) {
          cout << buffer2[i][j] << ",";
       }
       cout << endl;
    }
    */

    /*for (int i = 0; i < counter2; i ++) {
       err_pos_mag2 = 0.;
       for (int j = 0; j < 3; j ++) {
          err_pos_mag2 += error[i][j] * error[i][j];
       }
       err_pos_mag = sqrt(err_pos_mag2);
       cout << buffer2[i][9] << "," << err_pos_mag << endl; // (Time, magnitude of position error from interpolant)
    }
    */
    // Cannot delete G4TransportationManager::GetInstance();

    return 0;
}


  
