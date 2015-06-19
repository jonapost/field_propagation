// Nystrom stepper implemenation by Jason Suagee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 27 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license

#include "G4UniformMagField.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include <iomanip>


#include "G4QuadrupoleMagField.hh"
#include "G4CachedMagneticField.hh"

#include "G4CashKarpRKF45.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "BogackiShampine23.hh"
#include "G4LineSection.hh"
#include "G4MagIntegratorStepper.hh"
#include "DormandPrince745.hh"
#include "BogackiShampine45.hh"
#include "G4ClassicalRK4.hh"
#include "G4SimpleHeum.hh"
#include "G4ChargeState.hh"

#include "Mag_UsualEqRhs_IntegrateByTime.hh"

#include "ChawlaSharmaRKNstepper.hh"
#include "G4NystromRK4.hh"
#include "ChawlaSharmaRKNstepperFullRenormalization.hh"
#include "NystromRK4noNormalization.hh"
#include "ChawlaSharmaWrapper.hh"

using namespace std;
using namespace CLHEP;

//Version for testing of ChawlaSharmaRKNstepper (Jason S.).

/*
 Usage is ./modifiedtestH [stepper no.] [No. Steps] [Step len] [choice of Mag field]

 Choice of Mag field:
 1: Uniform Mag Field
 2: Quadropole Mag Field

 Stepper No.
 Not included in this version: 0: G4ExactHelixStepper
 1: G4CashKarpRKF45
 2: BogackiShampine23
 3: DormandPrince745
 4: BogackiShampine45
 5: G4ClassicalRK4
 6: G4SimpleHeum
 Not included in this version: 7: BogackiShampine23modified (No renormalization inside of a step.)
 8: ChawlaSharmaRKNstepper
 9: ChawlaSharmaRKNstepper (Renormalization in between steps.)
 10: ChawlaSharmaRKNstepperFullRenormalization
 (Normalization in between steps, and inside of steps.)
 11: G4NystromRK4
 12: NystromRK4noNormalization
 */

int main(int argc, char *args[]) {

   /* -----------------------------SETTINGS-------------------------------- */

   /* USER SPACE
    - Modify values here */

   int no_of_steps = 100;        //Default No. of Steps for the stepper
   int stepper_no = 3;           //Choose stepper no., for reference see above
   G4double step_len = 10 * mm;    //Step length in milimeters

   //Set coordinates here
   G4double x_pos = 0.,                   //pos - position
         y_pos = 0., z_pos = 0.,


         // Might want to play around with the momentum values:
         x_mom = 0.0,                   //mom - momentum
         y_mom = 0.6, z_mom = 0.8,

         x_field = 0.0 * tesla,         //Uniform Magnetic Field (x,y,z)
         y_field = 0.0 * tesla, z_field = -1. * tesla;



   G4double momentum = 0.5 * proton_mass_c2;
   G4double imom = 1. / sqrt( x_mom * x_mom + y_mom * y_mom + z_mom * z_mom );

   x_mom *= momentum * imom;
   y_mom *= momentum * imom;
   z_mom *= momentum * imom;


   G4double kineticEnergy =  momentum*momentum /
                     ( std::sqrt( momentum*momentum + proton_mass_c2 * proton_mass_c2 )
   		    + proton_mass_c2 );
          // G4double velocity = momentum / ( proton_mass_c2 + kineticEnergy );

   //Set Charge etc.
   G4double particleCharge = +1.0,  // in e+ units
         spin = 0.0,                        // ignore the spin
         magneticMoment = 0.0,             // ignore the magnetic moment
         mass = proton_mass_c2 + kineticEnergy;

   /*----------------------------END-SETTINGS-------------------------------*/

   /************************************XXXXXX*****************************************/

   /*-------------------------PREPARING STEPPER-----------------------------*/

   /* CODER SPACE
    - don't modify values here */

   //Checking for command line values :
   if (argc > 1)
      stepper_no = atoi(args[1]);
   if (argc > 2)
      no_of_steps = atoi(args[2]);
   if (argc > 3)
      step_len = (float) (atof(args[3]) * mm);

   //Initialising coordinates
   G4double yIn[] = { x_pos, y_pos, z_pos, x_mom/proton_mass_c2, y_mom/proton_mass_c2, z_mom/proton_mass_c2, 0., 0., 0., 0. };
   if (stepper_no != 8){
	   yIn[3] *= proton_mass_c2;
	   yIn[4] *= proton_mass_c2;
	   yIn[5] *= proton_mass_c2;
   }

   //Empty buckets for results
   G4double dydx[10]    =  { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
            yout[10]    =  { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
            yerr[10]    =  { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

   //1. Create a field :
   //Mag_UsualEqRhs_IntegrateByTime *fEquation;
   G4Mag_EqRhs *fEquation;
   G4UniformMagField *myUniformField;
   G4QuadrupoleMagField *quadrupoleMagField;
   G4CachedMagneticField *myQuadField;

   // Process args for type of Mag field:
   if (argc > 4) {
      G4int type_mag_field = atoi(args[4]);
      switch (type_mag_field) {
         case 1:
        	myUniformField = new G4UniformMagField(
                  G4ThreeVector(x_field, y_field, z_field));
        	if (stepper_no == 8)
        		fEquation = new Mag_UsualEqRhs_IntegrateByTime(myUniformField);
        	else
        		fEquation = new G4Mag_UsualEqRhs(myUniformField);
            break;
         case 2:
            quadrupoleMagField = new G4QuadrupoleMagField(
                  10. * tesla / (50. * cm));
            myQuadField = new G4CachedMagneticField(quadrupoleMagField,
                  1.0 * cm);
            fEquation = new Mag_UsualEqRhs_IntegrateByTime(myQuadField);
            break;
         default:
            cout
                  << "4th argument must be \"1\" for uniform or \"2\" for quadropole "
                  << endl;
            return 0;
      }
   } else {
      cout
            << "4th argument must be supplied and ust be \"1\" for uniform or \"2\" for quadropole "
            << endl;
      return 0;
   }

   //Create an Equation :
   //G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);

   G4ChargeState chargeState(particleCharge, // The charge can change (dynamic)
         spin = 0.0, magneticMoment = 0.0);

   fEquation->SetChargeMomentumMass(chargeState,
         G4ThreeVector(x_mom, y_mom, z_mom).mag(),
         //momentum magnitude
         mass);               //No place for mass in fEquation though


   //Create a stepper :
   G4MagIntegratorStepper * myStepper; // * exactStepper;

   G4SimpleHeum *simpleHeumStepper;
   ChawlaSharmaRKNstepper *myChawlaStepper;
   ChawlaSharmaRKNstepperFullRenormalization *myChawlaSharmaRKNstepperFullRenormalization;
   G4NystromRK4 *myNystromStepper;
   NystromRK4noNormalization *myNystromStepper_noNormalization;

   ChawlaSharmaWrapper *myCSWStepper;

   //Choose the stepper based on the command line argument
   switch (stepper_no) {
      // case 0:
      //	myStepper = new G4ExactHelixStepper(fEquation);
      //	break;
      case 1:
         myStepper = new G4CashKarpRKF45(fEquation);
         break;
      case 2:
         myStepper = new BogackiShampine23(fEquation);
         break;
      case 3:
         myStepper = new DormandPrince745(fEquation);
         break;
      case 4:
         myStepper = new BogackiShampine45(fEquation);
         break;
      case 5:
         myStepper = new G4ClassicalRK4(fEquation);
         break;
      case 6:
         myStepper = new G4SimpleHeum(fEquation);
         break;
      case 8:
         myStepper = new ChawlaSharmaRKNstepper(fEquation);
         break;
      case 9:
         myCSWStepper = new ChawlaSharmaWrapper(fEquation);
         break;
      case 10:
         myStepper = new ChawlaSharmaRKNstepperFullRenormalization(fEquation);
         // Normalization in between steps, and inside of steps.
         break;
      case 11:
         myStepper = new G4NystromRK4(fEquation);
         break;
      case 12:
         myStepper = new NystromRK4noNormalization(fEquation);
         break;
      default:
         myStepper = 0;
   }

   // For output to ipython notebook (for visualization)
   cout.setf(ios_base::fixed);
   cout.precision(10);


   // Probably could be done differently, but it's just a test program
   switch (stepper_no) {
      case 6:
         simpleHeumStepper = dynamic_cast<G4SimpleHeum*>(myStepper);
         break;
      case 8:
         //myChawlaStepper = dynamic_cast<ChawlaSharmaRKNstepper*>(myStepper);
         break;
      case 9:
         //myChawlaStepper = dynamic_cast<ChawlaSharmaRKNstepper*>(myStepper);
         break;
      case 10:
         myChawlaSharmaRKNstepperFullRenormalization =
               dynamic_cast<ChawlaSharmaRKNstepperFullRenormalization*>(myStepper);
         break;
      case 11:
         myNystromStepper = dynamic_cast<G4NystromRK4*>(myStepper);
         break;
      case 12:
         myNystromStepper_noNormalization =
               dynamic_cast<NystromRK4noNormalization*>(myStepper);
         break;
      default:
         break;
   }

   /*-----------------------END PREPARING STEPPER---------------------------*/

   /*----------------NOW STEPPING-----------------*/

   G4double mom_inorm;

   for (int j = 0; j < no_of_steps; j++) {
      switch (stepper_no) {
         case 6:
            simpleHeumStepper->RightHandSide(yIn, dydx);
            simpleHeumStepper->DumbStepper(yIn, dydx, step_len, yout);
            break;
         case 8:
            //myChawlaStepper->DumbStepper(yIn, step_len, yout);
        	myStepper->Stepper(yIn, dydx, step_len, yout, yerr);

        	 break;
         case 9:
            //myCSWStepper->Stepper(yIn, dydx, step_len, yout, yerr);
        	 myCSWStepper->DumbStepper(yIn, step_len, yout);
        	 break;
         case 10:
            // Renormalize momentum before the first step, and during each step.
            mom_inorm = 1.
                  / sqrt(yIn[3] * yIn[3] + yIn[4] * yIn[4] + yIn[5] * yIn[5]);
            for (int k = 3; k < 6; k++) {
               yIn[k] *= mom_inorm;
            }

            myChawlaSharmaRKNstepperFullRenormalization->DumbStepper(yIn,
                  step_len, yout);
            break;

         case 11:
            myNystromStepper->ComputeRightHandSide(yIn, dydx);
            myNystromStepper->Stepper(yIn, dydx, step_len, yout, yerr);
            break;
         case 12:
            // Normalize momentum before compute RHS
            mom_inorm = 1.
                  / sqrt(yIn[3] * yIn[3] + yIn[4] * yIn[4] + yIn[5] * yIn[5]);
            for (int k = 3; k < 6; k++) {
               yIn[k] *= mom_inorm;
            }

            myNystromStepper_noNormalization->ComputeRightHandSide(yIn, dydx);
            myNystromStepper_noNormalization->Stepper(yIn, dydx, step_len, yout,
                  yerr);
            break;
         default:
            myStepper->RightHandSide(yIn, dydx);
            myStepper->Stepper(yIn, dydx, step_len, yout, yerr); //call the stepper
            mom_inorm = 1. / sqrt(yIn[3] * yIn[3] + yIn[4] * yIn[4] + yIn[5] * yIn[5]);
			for (int k = 3; k < 6; k++) {
			   yIn[k] *= mom_inorm;
			}

      }

      // exactStepper->RightHandSide(yInX, dydxRef);
      // exactStepper->Stepper(yInX, dydxRef, step_len, youtX, yerrX);

      //-> Then print the data
      // Modified to just print out X,Y,Z position data,
      // so is easier to feed into ipython notebook.
      //for (int i = 0; i < 6; i++){
      for (int i = 0; i < 2; i++) {
         // Uncomment out if you want to print error
         // cout << yout[i] - youtX[i] << ",";
         cout << yout[i] << ",";
      }
      // Uncomment out if you want to print error
      //cout << yout[2] - youtX[2] << endl;
      cout << yout[2];

      //Copy yout into yIn
      for (int i = 0; i < 6; i++) {
         yIn[i] = yout[i];      }
      cout << "\n";
   }

   switch (stepper_no) {
      case 8:
         delete myStepper;
         break;
      case 9:
         delete myCSWStepper;
         break;
      case 10:
         delete myChawlaSharmaRKNstepperFullRenormalization;
         break;
      case 11:
         delete myNystromStepper;
         break;
      case 12:
         delete myNystromStepper_noNormalization;
         break;
      case 6:
         delete simpleHeumStepper;
         break;
      default:
         delete myStepper;
         break;

   }
   // delete exactStepper;
}
