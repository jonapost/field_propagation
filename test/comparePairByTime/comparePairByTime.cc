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

#include "MagIntegratorStepperbyTime.hh"
#include "FineRKNG34.hh"


#include <iostream>
#include "G4ThreeVector.hh"

using namespace std;
using namespace CLHEP;



int main(int argc, char *args[]) {
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

   //G4double imass = 1. / mass;

   G4double yIn1[10] = { x_pos, y_pos, z_pos, x_mom, y_mom, z_mom, 0., 0., 0., 0. };
   G4double dydx1[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

   G4double yIn2[10] = { x_pos, y_pos, z_pos, x_mom, y_mom, z_mom, 0., 0., 0., 0. };
   G4double dydx2[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

   G4UniformMagField *myUniformField;
   G4QuadrupoleMagField *quadrupoleMagField;
   G4Mag_EqRhs *fEquation;

   G4int mag_field_choice = 2;
   if (argc > 5)
      mag_field_choice = atoi(args[5]);

   if (mag_field_choice == 1){
      myUniformField = new G4UniformMagField(
                        G4ThreeVector(x_field, y_field, z_field));
      fEquation = new Mag_UsualEqRhs_IntegrateByTime(myUniformField);
   }
   else {
      quadrupoleMagField = new G4QuadrupoleMagField(
                        10. * tesla / (50. * cm));
      // myQuadField = new G4CachedMagneticField(quadrupoleMagField,
      //      1.0 * cm);
      fEquation = new Mag_UsualEqRhs_IntegrateByTime(quadrupoleMagField);
   }
   G4ChargeState chargeState(particleCharge, // The charge can change (dynamic)
               spin = 0.0, magneticMoment = 0.0);
   fEquation->SetChargeMomentumMass(chargeState,
         G4ThreeVector(x_mom, y_mom, z_mom).mag(),
         //momentum magnitude
         mass);

   G4int stepper_no1 = 0, stepper_no2 = 0, no_of_steps = 10;
   G4double step_len = 1.;


   if (argc > 1)
         stepper_no1 = atoi(args[1]);
   if (argc > 2)
            stepper_no2 = atoi(args[2]);
   if (argc > 3)
      no_of_steps = atoi(args[3]);
   if (argc > 4)
      step_len = (float) (atof(args[4]) * mm);

   G4MagIntegratorStepper *myStepper1;
   G4MagIntegratorStepper *myStepper2;

   switch (stepper_no1) {
      // case 0:
      // myStepper = new G4ExactHelixStepper(fEquation);
      // break;
      case 1:
         myStepper1 = new MagIntegratorStepper_byTime<G4CashKarpRKF45>(fEquation);
         break;
      case 2:
         myStepper1 = new MagIntegratorStepper_byTime<BogackiShampine23>(fEquation);
         break;
      case 3:
         myStepper1 = new MagIntegratorStepper_byTime<DormandPrince745>(fEquation);
         break;
      case 4:
         myStepper1 = new MagIntegratorStepper_byTime<BogackiShampine45>(fEquation);
         break;
      case 5:
         myStepper1 = new MagIntegratorStepper_byTime<G4ClassicalRK4>(fEquation);
         break;
      case 6:
         myStepper1 = new MagIntegratorStepper_byTime<G4SimpleHeum>(fEquation);
         break;
      case 7:
         myStepper1 = new MagIntegratorStepper_byTime<G4NystromRK4>(fEquation);
         break;
      case 8:
         myStepper1 = new MagIntegratorStepper_byTime<ChawlaSharmaRKNstepper>(fEquation);
         break;
      case 9:
         myStepper1 = new MagIntegratorStepper_byTime<FineRKNG34>(fEquation);
         break;
      default:
         myStepper1 = 0;
   }
   switch (stepper_no2) {
      // case 0:
      // myStepper = new G4ExactHelixStepper(fEquation);
      // break;
      case 1:
         myStepper2 = new MagIntegratorStepper_byTime<G4CashKarpRKF45>(fEquation);
         break;
      case 2:
         myStepper2 = new MagIntegratorStepper_byTime<BogackiShampine23>(fEquation);
         break;
      case 3:
         myStepper2 = new MagIntegratorStepper_byTime<DormandPrince745>(fEquation);
         break;
      case 4:
         myStepper2 = new MagIntegratorStepper_byTime<BogackiShampine45>(fEquation);
         break;
      case 5:
         myStepper2 = new MagIntegratorStepper_byTime<G4ClassicalRK4>(fEquation);
         break;
      case 6:
         myStepper2 = new MagIntegratorStepper_byTime<G4SimpleHeum>(fEquation);
         break;
      case 7:
         myStepper2 = new MagIntegratorStepper_byTime<G4NystromRK4>(fEquation);
         break;
      case 8:
         myStepper2 = new MagIntegratorStepper_byTime<ChawlaSharmaRKNstepper>(fEquation);
         break;
      case 9:
         myStepper2 = new MagIntegratorStepper_byTime<FineRKNG34>(fEquation);
         break;
      default:
         myStepper2 = 0;
   }
   /*
   if ( ! myStepper->mass ){
      cout << "no mass" << endl;
      return 0;
   }
   */
   // MagIntegratorStepper_byTime<G4ClassicalRK4>  *myStepper = new MagIntegratorStepper_byTime<G4ClassicalRK4>(fEquation);

   // For output to ipython notebook (for visualization)
   cout.setf(ios_base::fixed);
   cout.precision(25);

   /*-----------------------END PREPARING STEPPER---------------------------*/

   /*----------------NOW STEPPING-----------------*/

   G4double yout1[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   G4double yerr1[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

   G4double yout2[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   G4double yerr2[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };



   /*
   for (int k = 0; k < 3; k++) {
      cout << yIn[k] << ",";
   }
   for (int k = 3; k < 6; k++) {
      cout << yIn[k] / mass << ",";
   }
   cout << endl;
   cout << endl;
   */

   G4double error;

   for (int j = 0; j < no_of_steps; j++) {
      myStepper1->ComputeRightHandSide(yIn1, dydx1);
      myStepper1->Stepper(yIn1, dydx1, step_len, yout1, yerr1); //call the stepper
      myStepper2->ComputeRightHandSide(yIn2, dydx2);
      myStepper2->Stepper(yIn2, dydx2, step_len, yout2, yerr2); //call the stepper

      // Error Magnitude output:
      error = 0.;
      for (int k = 0; k < 3; k++)
         error += ( yout1[k] - yout2[k] ) * ( yout1[k] - yout2[k] );
      cout << sqrt(error);

      // Position output:
      /* for (int k = 0; k < 3; k++) {
         cout << yout[k] << ",";
      }
      // Velocity output
      for (int k = 3; k < 6; k++) {
         // Uncomment to print out momentums:
         cout << yout[k] / mass << ",";
      }
      */
      cout << endl;

      //Copy yout into yIn
      for (int k = 0; k < 8; k++){
         yIn1[k] = yout1[k];
         yIn2[k] = yout2[k];

      }
   }

   delete myStepper1;
   delete myStepper2;
   delete fEquation;

}




