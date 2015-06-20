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

#include "MagIntegratorStepperbyTime.hh"

using namespace std;
using namespace CLHEP;


MagIntegratorStepper_byTime<G4CashKarpRKF45>  *myCashKarp;
MagIntegratorStepper_byTime<BogackiShampine23> *myBogacki23;
MagIntegratorStepper_byTime<DormandPrince745> *myDormand745;
MagIntegratorStepper_byTime<BogackiShampine45> *myBogacki45;
MagIntegratorStepper_byTime<G4ClassicalRK4> *myClassicRK;
MagIntegratorStepper_byTime<G4SimpleHeum> *mySimpleHeum;
MagIntegratorStepper_byTime<G4NystromRK4> *myNystrom;
MagIntegratorStepper_byTime<ChawlaSharmaRKNstepper> *myCharma;


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


   G4double yIn[10] = { x_pos, y_pos, z_pos, x_mom, y_mom, z_mom, 0., 0., 0., 0. };
   G4double dydx[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };


   G4UniformMagField *myUniformField = new G4UniformMagField(
                     G4ThreeVector(x_field, y_field, z_field));
   G4Mag_EqRhs *fEquation = new G4Mag_UsualEqRhs(myUniformField);

   G4ChargeState chargeState(particleCharge, // The charge can change (dynamic)
               spin = 0.0, magneticMoment = 0.0);
   fEquation->SetChargeMomentumMass(chargeState,
         G4ThreeVector(x_mom, y_mom, z_mom).mag(),
         //momentum magnitude
         mass);

   G4int stepper_no = 0, no_of_steps = 10;
   G4double step_len = 1.;


   if (argc > 1)
         stepper_no = atoi(args[1]);
   if (argc > 2)
      no_of_steps = atoi(args[1]);
   if (argc > 3)
      step_len = (float) (atof(args[2]) * mm);


   switch (stepper_no) {
            // case 0:
            // myStepper = new G4ExactHelixStepper(fEquation);
            // break;
            case 1:
               myCashKarp = new MagIntegratorStepper_byTime<G4CashKarpRKF45>(fEquation);
               break;
            case 2:
               myBogacki23 = new MagIntegratorStepper_byTime<BogackiShampine23>(fEquation);
               break;
            case 3:
               myDormand745 = new MagIntegratorStepper_byTime<DormandPrince745>(fEquation);
               break;
            case 4:
               myBogacki45 = new MagIntegratorStepper_byTime<BogackiShampine45>(fEquation);
               break;
            case 5:
               myClassicRK = new MagIntegratorStepper_byTime<G4ClassicalRK4>(fEquation);
               break;
            case 6:
               mySimpleHeum = new MagIntegratorStepper_byTime<G4SimpleHeum>(fEquation);
               break;
            case 7:
               myNystrom = new MagIntegratorStepper_byTime<G4NystromRK4>(fEquation);
               break;
            case 8:
               myCharma = new MagIntegratorStepper_byTime<ChawlaSharmaRKNstepper>(fEquation);
               break;
            default:
               myStepper = 0;
         }
      }





   MagIntegratorStepper_byTime<G4ClassicalRK4>  *myStepper = new MagIntegratorStepper_byTime<G4ClassicalRK4>(fEquation);

   // For output to ipython notebook (for visualization)
   cout.setf(ios_base::fixed);
   cout.precision(10);

   /*-----------------------END PREPARING STEPPER---------------------------*/

   /*----------------NOW STEPPING-----------------*/

   G4double yout[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   G4double yerr[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };


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

   for (int j = 0; j < no_of_steps; j++) {
      myStepper->ComputeRightHandSide(yIn, dydx);
      myStepper->Stepper(yIn, dydx, step_len, yout, yerr); //call the stepper


      // Position output:
      for (int k = 0; k < 3; k++) {
         cout << yout[k] << ",";
      }
      // Velocity output
      for (int k = 3; k < 6; k++) {
         // Uncomment to print out momentums:
         cout << yout[k] / mass << ",";
      }
      cout << endl;

      //Copy yout into yIn
      for (int k = 0; k < 8; k++){
         yIn[k] = yout[k];
      }
   }


   delete myStepper;
   delete fEquation;

}




