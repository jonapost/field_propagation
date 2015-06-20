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

#include "G4NystromRK4.hh"

#include "Mag_UsualEqRhs_IntegrateByTime.hh"
#include "ChawlaSharmaRKNstepper.hh"
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
 7: G4NystromRK4
 8: ChawlaSharmaWrapper
 */

int main(int argc, char *args[]) {

   /* -----------------------------SETTINGS-------------------------------- */

   /* USER SPACE
    - Modify values here */

   int no_of_steps = 100;        //Default No. of Steps for the stepper
   int stepper_choices[2];       //Choose stepper no., for reference see above
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
      stepper_choices[0] = atoi(args[1]);
   if (argc > 2)
      stepper_choices[1] = atoi(args[2]);
   if (argc > 3)
      no_of_steps = atoi(args[3]);
   if (argc > 4)
      step_len = (float) (atof(args[4]) * mm);

   //Initialising coordinates
   G4double **yIn = new G4double*[2];

   G4double Begining_vals[10] = { x_pos, y_pos, z_pos, x_mom, y_mom, z_mom, 0., 0., 0., 0. };
   G4double zeros[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };


   for (int i = 0; i < 2; i ++){
      yIn[i] = new G4double[10];
      for (int j = 0; j < 10; j ++)
         yIn[i][j] = Begining_vals[j];
   }

   //Empty buckets for results
   G4double **dydx = new G4double*[2];
   G4double **yout = new G4double*[2];
   G4double **yerr = new G4double*[2];

   for (int i = 0; i < 2; i ++){
      dydx[i] = new G4double[10];
      yout[i] = new G4double[10];
      yerr[i] = new G4double[10];
      for (int j = 0; j < 10; j ++) {
         dydx[i][j] = zeros[j];
         yout[i][j] = zeros[j];
         yerr[i][j] = zeros[j];
      }
   }
   //1. Create a field :
   //Mag_UsualEqRhs_IntegrateByTime *fEquation;
   G4Mag_EqRhs **fEquation = new G4Mag_EqRhs*[2];
   G4UniformMagField *myUniformField;
   G4QuadrupoleMagField *quadrupoleMagField;
   G4CachedMagneticField *myQuadField;

   // Process args for type of Mag field:
   G4int type_mag_field;
   if (argc > 5)
      type_mag_field = atoi(args[5]);
   if (type_mag_field) {
      switch (type_mag_field) {
         case 1:
        	myUniformField = new G4UniformMagField(
                  G4ThreeVector(x_field, y_field, z_field));
        	if (stepper_choices[0] == 8)
        		fEquation[0] = new Mag_UsualEqRhs_IntegrateByTime(myUniformField);
        	else
        		fEquation[0] = new G4Mag_UsualEqRhs(myUniformField);
         if (stepper_choices[1] == 8)
            fEquation[1] = new Mag_UsualEqRhs_IntegrateByTime(myUniformField);
         else
            fEquation[1] = new G4Mag_UsualEqRhs(myUniformField);
         break;

         case 2:
            quadrupoleMagField = new G4QuadrupoleMagField(
                  10. * tesla / (50. * cm));
            myQuadField = new G4CachedMagneticField(quadrupoleMagField,
                  1.0 * cm);
            if (stepper_choices[0] == 8)
               fEquation[0] = new Mag_UsualEqRhs_IntegrateByTime(myQuadField);
            else
               fEquation[0] = new G4Mag_UsualEqRhs(myQuadField);
            if (stepper_choices[1] == 8)
               fEquation[1] = new Mag_UsualEqRhs_IntegrateByTime(myQuadField);
            else
               fEquation[1] = new G4Mag_UsualEqRhs(myQuadField);
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

   G4ChargeState chargeState1(particleCharge, // The charge can change (dynamic)
         spin = 0.0, magneticMoment = 0.0);

   G4ChargeState chargeState2(particleCharge, // The charge can change (dynamic)
            spin = 0.0, magneticMoment = 0.0);


   fEquation[0]->SetChargeMomentumMass(chargeState1,
         G4ThreeVector(x_mom, y_mom, z_mom).mag(),
         //momentum magnitude
         mass);               //No place for mass in fEquation though
   fEquation[1]->SetChargeMomentumMass(chargeState2,
            G4ThreeVector(x_mom, y_mom, z_mom).mag(),
            //momentum magnitude
            mass);               //No place for mass in fEquation though


   //Create a stepper :
   G4MagIntegratorStepper **myStepper = new G4MagIntegratorStepper*[2]; // a length 2 array to hold our steppers we want to compare

   //G4SimpleHeum *simpleHeumStepper;
   G4NystromRK4 *myNystromStepper;
   ChawlaSharmaWrapper *myCSWStepper;

   //Choose the stepper based on the command line argument

   G4int stepper_no;
   for (int i = 0; i < 2; i ++) {
      stepper_no = stepper_choices[i];

      switch (stepper_no) {
         // case 0:
         //	myStepper = new G4ExactHelixStepper(fEquation);
         //	break;
         case 0:
            if (type_mag_field != 1){
               cout << "must use Uniform Mag field with Exact Helix Stepper!" << endl;
               return 0;
            }
            else {
               myStepper[i] = new G4ExactHelixStepper(fEquation[i]);
               break;
            }
         case 1:
            myStepper[i] = new G4CashKarpRKF45(fEquation[i]);
            break;
         case 2:
            myStepper[i] = new BogackiShampine23(fEquation[i]);
            break;
         case 3:
            myStepper[i] = new DormandPrince745(fEquation[i]);
            break;
         case 4:
            myStepper[i] = new BogackiShampine45(fEquation[i]);
            break;
         case 5:
            myStepper[i] = new G4ClassicalRK4(fEquation[i]);
            break;
         case 6:
            myStepper[i] = new G4SimpleHeum(fEquation[i]);
            break;
         case 7:
            myNystromStepper = new G4NystromRK4(fEquation[i]);
            break;
         case 8:
            myCSWStepper = new ChawlaSharmaWrapper(fEquation[i]);
            break;
         default:
            myStepper[i] = 0;
      }
   }


   // For output to ipython notebook (for visualization)
   cout.setf(ios_base::fixed);
   cout.precision(10);

   /*-----------------------END PREPARING STEPPER---------------------------*/

   /*----------------NOW STEPPING-----------------*/

   G4double total_time[2] = {0., 0.};
   G4double total_distance[2] = {0., 0.};
   G4double stored_vals[10] = { 0., 0., 0., 0., 0., 0., 0., 0.};
   G4double difference_values[10] = { 0., 0., 0., 0., 0., 0., 0., 0.};
   G4double mag = 0.;

   G4String output_format;
   if (argc > 6)
      output_format = args[6];

   for (int j = 0; j < no_of_steps; j++) {
      for (int i = 0; i < 2; i ++){
         stepper_no = stepper_choices[i];
         switch (stepper_no) {
            /*case 6:
               simpleHeumStepper->RightHandSide(yIn, dydx);
               simpleHeumStepper->DumbStepper(yIn, dydx, step_len, yout);
               break;
            */
            case 7:
               myNystromStepper->ComputeRightHandSide(yIn[i], dydx[i]);
               myNystromStepper->Stepper(yIn[i], dydx[i], step_len, yout[i], yerr[i]);
               break;
            case 8:
               //myCSWStepper->Stepper(yIn, dydx, step_len, yout, yerr);
               myCSWStepper->DumbStepper(yIn[i], step_len, yout[i]);
               break;

            default:
               myStepper[i]->RightHandSide(yIn[i], dydx[i]);
               myStepper[i]->Stepper(yIn[i], dydx[i], step_len, yout[i], yerr[i]); //call the stepper
               //mom_inorm = 1. / sqrt(yIn[3] * yIn[3] + yIn[4] * yIn[4] + yIn[5] * yIn[5]);
         }

         if (output_format == G4String("values")) {
            // Position output:
            for (int k = 0; k < 3; k++) {
               cout << yout[i][k] << ",";
            }

            // Velocity output
            for (int k = 3; k < 6; k++) {
               // Uncomment to print out momentums:
               //cout << yout[i][k] / mass << ",";
            }


            // Total time:
            if (stepper_no == 8)
               total_time[i] = yout[i][6];
            else
               total_time[i] += step_len / ( G4ThreeVector( yIn[i][3], yIn[i][4], yIn[i][5] ).mag() / mass ); // Step / ||v||. Remember that yIn[3..5] represents momentum.
            // Uncomment to print out time
            //cout << total_time[i] << ",";


            // Total Arc Length:
            if (stepper_no == 8)
               total_distance[i] = yout[i][7];
            else
               total_distance[i] += step_len;
            // Uncomment to print arc_length
            //cout << total_distance[i];

            if (i == 0)
               cout << ",";
            else
               cout << endl;
         }
         else {
            if (i == 0) {
               for (int k = 0; k < 8; k ++)
                  stored_vals[k] = yout[i][k];
            }
            else {
               for (int k = 0; k < 8; k ++)
                  difference_values[k] = yout[i][k] - yout[i-1][k]; // stored_vals[k];

               mag = 0.;
               for (int k = 0; k < 3; k ++)
                  mag += difference_values[k] * difference_values[k];
               cout << sqrt(mag) << ",";

               mag = 0.;
               for (int k = 3; k < 6; k ++)
                  mag += difference_values[k] * difference_values[k];
               cout << sqrt(mag);
               cout << endl;
            }
         }
         //Copy yout into yIn
         for (int k = 0; k < 8; k++){
            yIn[i][k] = yout[i][k];
         }
      }

   }

   for (int i = 0; i < 2; i ++) {
      switch (stepper_choices[i]) {
         case 8:
            delete myCSWStepper;
            break;
         case 7:
            delete myNystromStepper;
            break;
         default:
            delete myStepper[i];
            break;
      }
      delete fEquation[i];
   }
   delete[] dydx;
   delete[] yIn;
   delete[] yout;
   delete[] yerr;
   delete[] myStepper;
   delete[] fEquation;
}

