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


#include "G4MagHelicalStepper.hh"

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
#include "MagEqRhsbyTimestoreB.hh"


#include "ChawlaSharmaRKNstepper.hh"

#include "MagIntegratorStepperbyTime.hh"
#include "MagIntegratorStepperByArcLength.hh"

#include "FineRKNG34.hh"
#include "FineRKNG45.hh"
#include "MuruaRKN6459.hh"
#include "MuruaRKN5459.hh"

//#include "VernerRK78.hh"
//#include "TsitourasRK45.hh"

#include <iostream>
#include "G4ThreeVector.hh"


#include "Interpolant.hh"

#include "isTracking.hh"

#ifdef TRACKING
#include "StepTracker.hh"
#endif

using namespace std;
using namespace CLHEP;



int main(int argc, char *args[]) {

   //Position = G4ThreeVector(-100.*mm, 50.*mm, 150.*mm)


   G4double x_pos = 0.,                   //pos - position
            y_pos = 0., z_pos = 0.,
   //G4double x_pos = -100.*mm,                   //pos - position
   //         y_pos = 50.*mm, z_pos = 150.*mm,




            // Might want to play around with the momentum values:
            x_mom = 0.2,                   //mom - momentum
            y_mom = 0.6, z_mom = 0.8,

            x_field = 0.0 * tesla,         //Uniform Magnetic Field (x,y,z)
            y_field = 0.0 * tesla, z_field = -1. * tesla;


   G4double momentum = 0.5 * proton_mass_c2;
   G4double imom = 1. / sqrt( x_mom * x_mom + y_mom * y_mom + z_mom * z_mom );

   //x_mom *= momentum * imom;
   //y_mom *= momentum * imom;
   //z_mom *= momentum * imom;

   x_mom *= imom;
   y_mom *= imom;
   z_mom *= imom;


   G4double kineticEnergy =  momentum*momentum /
                           ( std::sqrt( momentum*momentum + proton_mass_c2 * proton_mass_c2 )
                             + proton_mass_c2 );
          // G4double velocity = momentum / ( proton_mass_c2 + kineticEnergy );

   //Set Charge etc.
   G4double particleCharge = +1.0,  // in e+ units
         spin = 0.0,                        // ignore the spin
         magneticMoment = 0.0,             // ignore the magnetic moment
         mass = proton_mass_c2; // + kineticEnergy;

   //G4double imass = 1. / mass;



   G4double coeff = std::sqrt( kineticEnergy*kineticEnergy + 2.*mass*kineticEnergy );



   G4double yIn[10] = { x_pos, y_pos, z_pos, coeff*x_mom, coeff*y_mom, coeff*z_mom, 0., 0., 0., 0. };
   G4double dydx[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

   //G4UniformMagField *myUniformField;
   //G4QuadrupoleMagField *quadrupoleMagField;
   G4MagneticField *my_non_cached_MagField;
   G4CachedMagneticField *myMagField;
   G4Mag_EqRhs *fEquation;

   G4int mag_field_choice = 2; // Choice is quadropole mag field
   // if (argc > 4)
   //   mag_field_choice = atoi(args[4]);

   if (mag_field_choice == 1){
      my_non_cached_MagField = new G4UniformMagField(
                        G4ThreeVector(x_field, y_field, z_field));
      myMagField = new G4CachedMagneticField( my_non_cached_MagField, 0.0 * cm);

      //fEquation = new Mag_UsualEqRhs_IntegrateByTime(myUniformField);
   }
   else {
      my_non_cached_MagField = new G4QuadrupoleMagField(
                        10. * tesla / (50. * cm));
      myMagField = new G4CachedMagneticField( my_non_cached_MagField, 0.0 * cm);
      // myQuadField = new G4CachedMagneticField(quadrupoleMagField,
      //      1.0 * cm);
      //fEquation = new Mag_UsualEqRhs_IntegrateByTime(quadrupoleMagField);
   }


   G4int stepper_no = 1, no_of_steps = 30;
   G4double step_len = 10.;
   char *outfile_name, *meta_outfile_name;

   if (argc > 1)
      stepper_no = atoi(args[1]);

   if (argc > 2)
      step_len = (float) (atof(args[2]) * mm);

   if (argc > 3)
      no_of_steps = atoi(args[3]);

   if (argc > 4)
      outfile_name = args[4];
   else
      outfile_name = "out";

   if (argc > 5)
      meta_outfile_name= args[5];
   else
      meta_outfile_name = "meta_out";

   G4MagIntegratorStepper *pStepper;
   //G4MagHelicalStepper *pHelicalStepper;

      //G4cout << " Setting up field of type: " << fieldName << G4endl;
   switch ( stepper_no )
   {
      /*
      case -3: // MuruaRKN6459
         //fEquation = new MagEqRhs_byTime_storeB(myMagField);

         fEquation = new G4Mag_UsualEqRhs(myMagField);
         pStepper = new MagIntegratorStepper_byArcLength<MuruaRKN6459>( fEquation );
         break;
      case -2:
         fEquation = new G4Mag_UsualEqRhs(myMagField);
         pStepper = new MagIntegratorStepper_byArcLength<FineRKNG45>( fEquation );
         break;
      */
      //case -1:
      //   fEquation = new G4Mag_UsualEqRhs(myMagField);
      //   pStepper = new MagIntegratorStepper_byArcLength<FineRKNG34>( fEquation );
      //   break;
      case -1:
         //fEquation = new MagEqRhs_byTime_storeB(myMagField);

         fEquation = new MagEqRhs_byTime_storeB(myMagField);
         pStepper = new MagIntegratorStepper_byTime<MuruaRKN5459>( fEquation );
         break;

      case 0: // MuruaRKN6459
         //fEquation = new MagEqRhs_byTime_storeB(myMagField);

         fEquation = new Mag_UsualEqRhs_IntegrateByTime(myMagField);
         pStepper = new MagIntegratorStepper_byTime<MuruaRKN6459>( fEquation );
         break;
      case 1:
         fEquation = new Mag_UsualEqRhs_IntegrateByTime(myMagField);
         pStepper = new MagIntegratorStepper_byTime<FineRKNG45>( fEquation );
         break;
      case 2:
         fEquation = new Mag_UsualEqRhs_IntegrateByTime(myMagField);
         pStepper = new MagIntegratorStepper_byTime<FineRKNG34>( fEquation );
         break;
      case 3:
         fEquation = new Mag_UsualEqRhs_IntegrateByTime(myMagField);
         pStepper = new MagIntegratorStepper_byTime<ChawlaSharmaRKNstepper>( fEquation );
         break;
      case 4:
         fEquation = new G4Mag_UsualEqRhs(myMagField);
         pStepper = new MagIntegratorStepper_byArcLength<G4CashKarpRKF45>( fEquation );
         break;
      case 5:
         fEquation = new G4Mag_UsualEqRhs(myMagField);
         pStepper = new MagIntegratorStepper_byArcLength<G4ClassicalRK4>( fEquation );
         break;
      case 6:
         fEquation = new G4Mag_UsualEqRhs(myMagField);
         pStepper = new MagIntegratorStepper_byArcLength<G4SimpleHeum>( fEquation );
         break;
      case 7:
         fEquation = new G4Mag_UsualEqRhs(myMagField);
         pStepper = new MagIntegratorStepper_byArcLength<DormandPrince745>( fEquation );
         break;

      case 8:
         fEquation = new G4Mag_UsualEqRhs(myMagField);
         pStepper = new MagIntegratorStepper_byArcLength<BogackiShampine45>( fEquation );
         break;


      //case 9:
      //   fEquation = new G4Mag_UsualEqRhs(myMagField);
      //   pHelicalStepper = new G4ExactHelixStepper(fEquation);

   }


   G4ChargeState chargeState(particleCharge, // The charge can change (dynamic)
                             spin = 0.0, magneticMoment = 0.0);
   fEquation->SetChargeMomentumMass(chargeState,
                                    G4ThreeVector(x_mom, y_mom, z_mom).mag(),
                                    //momentum magnitude
                                    mass);


#ifdef TRACKING

   pStepper -> ComputeRightHandSide(yIn, dydx);
   G4double beginning[11] = { 0., 0., x_pos, y_pos, z_pos, x_mom/mass, y_mom/mass, z_mom/mass, dydx[3]/mass, dydx[4]/mass, dydx[5]/mass };

   //G4double beginning[11] = { 0., 0., x_pos, y_pos, z_pos, x_mom/mass, y_mom/mass, z_mom/mass, dydx[3], dydx[4], dydx[5] };

   StepTracker *myStepTracker = new StepTracker( beginning );

   pStepper -> setTracker(myStepTracker);


   myStepTracker -> set_within_AdvanceChordLimited(true);

#endif

   // MagIntegratorStepper_byTime<G4ClassicalRK4>  *myStepper = new MagIntegratorStepper_byTime<G4ClassicalRK4>(fEquation);

   // For output to ipython notebook (for visualization)
   cout.setf(ios_base::fixed);
   cout.precision(15);

   /*-----------------------END PREPARING STEPPER---------------------------*/

   /*----------------NOW STEPPING-----------------*/

   G4double yout[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   G4double yerr[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

   //G4double yout_copy[10] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

   //G4double y_interpolated[3];

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


   //Interpolant *minterpolant = new Interpolant();

   //G4double error;

   for (int j = 0; j < no_of_steps; j++) {


      //if (! (stepper_no == 9) ) {

      pStepper->ComputeRightHandSide(yIn, dydx);

#ifdef TRACKING
      myStepTracker -> ArmTracker();
#endif

      pStepper->Stepper(yIn, dydx, step_len, yout, yerr); //call the stepper

      //}

      /*
      if (stepper_no == 9) {

         pHelicalStepper -> AdvanceHelix(yIn, G4ThreeVector(x_field, y_field, z_field), step_len, yout);

         //pHelicalStepper -> ComputeRightHandSide(yIn, dydx);
         //pHelicalStepper -> Stepper(yIn, dydx, step_len, yout, yerr);

         //cout << pHelicalStepper -> GetRadHelix() << endl;
         return 0;

      }
      */





#ifdef TRACKING
      myStepTracker -> update_time_arclength(step_len / ( myStepTracker -> last_velocity() ), step_len );
#endif
      //Copy yout into yIn
      for (int k = 0; k < 8; k++){
         yIn[k] = yout[k];
      }

      // Output yerr:
      //for (int k = 0; k < 6; k ++)
      //   cout << yout[k] << ", ";
      //cout << endl;

   }
   // Record final step:

   /* Don't need this anymore:
   pStepper -> ComputeRightHandSide(yout, dydx);

   for (int i = 3; i < 6; i ++)
               yout[i] /= fEquation -> FMass();
   for (int i = 3; i < 6; i ++)
                  dydx[i] /= fEquation -> FMass();


   myStepTracker -> RecordResultOfStepper(yout, dydx);

   myStepTracker -> update_time_arclength(step_len / ( myStepTracker -> last_velocity() ), step_len );
   */

   //Output:

#ifdef TRACKING
   myStepTracker -> outputBuffer(outfile_name, meta_outfile_name);
#endif
   // Cleanup

   delete pStepper;

#ifdef TRACKING
   delete myStepTracker;
#endif
   delete fEquation;

   return 0;
}




