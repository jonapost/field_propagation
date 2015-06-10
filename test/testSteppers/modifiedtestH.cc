#include "G4UniformMagField.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4CashKarpRKF45.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "BogackiShampine23.hh"
#include "G4LineSection.hh"
#include "G4MagIntegratorStepper.hh"
#include "DormandPrince745.hh"
#include "BogackiShampine45.hh"
#include <iomanip>
#include "G4ClassicalRK4.hh"
#include "G4SimpleHeum.hh"
#include "G4ChargeState.hh"

#include "ChawlaSharmaRKNstepper.hh"
#include "BogackiShampine23modified.hh"
#include "G4NystromRK4.hh"
#include "ChawlaSharmaRKNstepperFullRenormalization.hh"


using namespace std;
using namespace CLHEP;

//Version for testing of ChawlaSharmaRKNstepper (Jason S.).

/* Stepper No.
   0: G4ExactHelixStepper
   1: G4CashKarpRKF45
   2: BogackiShampine23
   3: DormandPrince745
   4: BogackiShampine45
   5: G4ClassicalRK4
   6: G4SimpleHeum
   7: BogackiShampine23modified (No renormalization inside of a step.)
   8: ChawlaSharmaRKNstepper
   9: ChawlaSharmaRKNstepper (Renormalization in between steps.)
   10: ChawlaSharmaRKNstepperFullRenormalization 
       (Normalization in between steps, and inside of steps.)
   11: G4NystromRK4
*/

int main(int argc, char *args[]){

/* -----------------------------SETTINGS-------------------------------- */

/* USER SPACE
    - Modify values here */

   int no_of_steps = 100;        //Default No. of Steps for the stepper
   int stepper_no = 3;           //Choose stepper no., for refernce see above
   G4double step_len = 10*mm;    //Step length in milimeters
   
   //Set coordinates here
   G4double
   x_pos = 0.,                   //pos - position
   y_pos = 0.,
   z_pos = 0.,
   
   x_mom = 4.,                   //mom - momentum
   y_mom = 10.,
   z_mom = 2.,
   
   x_field = 0.05*tesla,         //Uniform Magnetic Field (x,y,z)
   y_field = 0.2*tesla,
   z_field = -0.1*tesla ;
   
   //Set Charge etc.
   G4double particleCharge = +1.0,  // in e+ units
   spin=0.0,                        // ignore the spin
   magneticMoment= 0.0,             // ignore the magnetic moment
   mass = 1;
   
   //Choice of output coordinates
   int
   columns[] =
      {
         1 ,  //x_pos
         1 ,  //y_pos
         1 ,  //z_pos
         0 ,  //x_mom
         0 ,  //y_mom
         0    //z_mom
      }; //Variables in yOut[] you want to display - 0 for No, 1 for yes
    
    /*----------------------------END-SETTINGS-------------------------------*/

/************************************XXXXXX*****************************************/

    
    
    
    /*-------------------------PREPARING STEPPER-----------------------------*/
    
/* CODER SPACE
    - don't modify values here */
    
    
   //Checking for command line values :
   if(argc>1)
      stepper_no = atoi(args[1]);
   if(argc > 2)
      no_of_steps = atoi(args[2]);
   if(argc > 3)
      step_len = (float)(atof(args[3])*mm);
   //Initialising coordinates
   G4double yIn[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
   
   //Empty buckets for results
   G4double dydx[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
   yout[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
   yerr[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   
   //1. Create a field :
   G4UniformMagField myField(G4ThreeVector(x_field, y_field, z_field));
   
   //Create an Equation :
   G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);
   
   G4ChargeState chargeState( particleCharge,      // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0);
   
   fEquation->SetChargeMomentumMass(chargeState,
                              G4ThreeVector(x_field, y_field, z_field).mag(), 
                                                   //momentum magnitude
                              mass);               //No place fo mass in fEquation though
   
   //Create a stepper :
   G4MagIntegratorStepper *myStepper; // , *exactStepper;
   ChawlaSharmaRKNstepper *myChawlaStepper;
   BogackiShampine23modified *myBogackiShampine23modified;
   ChawlaSharmaRKNstepperFullRenormalization 
                             *myChawlaSharmaRKNstepperFullRenormalization;
   G4NystromRK4 *myNystromStepper;
   
   //Choose the stepper based on the command line argument
   switch(stepper_no){
      case 0:
         myStepper = new G4ExactHelixStepper(fEquation);
         break;
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
      case 7:
         myStepper = new BogackiShampine23modified(fEquation);
         // No renormalization inside of a step.
         break;
      case 8:
         myStepper = new ChawlaSharmaRKNstepper(fEquation);
         break;
      case 9:
         myStepper = new ChawlaSharmaRKNstepper(fEquation);
         // Renormalization in between steps.
         break;
      case 10: 
         myChawlaSharmaRKNstepperFullRenormalization =
                     new ChawlaSharmaRKNstepperFullRenormalization(fEquation);
         // Normalization inbetween steps, and inside of steps.
         break;
      case 11:
         myStepper = new G4NystromRK4(fEquation);
         break;
      default :
         myStepper = 0 ;
    }
    
    // For output to Octave (for visualization)
    cout.setf (ios_base::fixed);
    cout.precision(9);
    
/*-----------------------END PREPARING STEPPER---------------------------*/
    

    
/*----------------NOW STEPPING-----------------*/
    
    
   switch(stepper_no){
      case 7:
         myBogackiShampine23modified = dynamic_cast<BogackiShampine23modified*>(myStepper);
         break;
      case 8:
         myChawlaStepper = dynamic_cast<ChawlaSharmaRKNstepper*>(myStepper);
         break;
      case 9:
         myChawlaStepper = dynamic_cast<ChawlaSharmaRKNstepper*>(myStepper);
         break;
      case 10:
         myChawlaSharmaRKNstepperFullRenormalization = 
            dynamic_cast<ChawlaSharmaRKNstepperFullRenormalization*>(myStepper);
         break;
      case 11:
         myNystromStepper = dynamic_cast<G4NystromRK4*>(myStepper);
         break;
      default:
         break;
   }

   G4double mom_norm;

   for(int j=0; j<no_of_steps; j++){
      switch(stepper_no){
         case 7:
            //myBogackiShampine23modified->modifiedRightHandSide(yIn, dydx);
            // Temp change to test a theory
            
            myBogackiShampine23modified->RightHandSide(yIn, dydx); 
            // Renormalize between steps ( Standard RightHandSide() )
            myBogackiShampine23modified->Stepper(yIn, dydx, step_len, yout, yerr);
            break;
         case 8:
            myChawlaStepper->DumbStepper(yIn, step_len, yout);
            break;
         case 9:
            // Renormalize momentum before the first step.
            mom_norm = sqrt(yIn[3]*yIn[3] + yIn[4]*yIn[4] + yIn[5]*yIn[5]);
            for (int k = 3; k < 6; k ++){
               yIn[k] /= mom_norm;
            }
            myChawlaStepper->DumbStepper(yIn, step_len, yout);
            break;
         case 10:
            // Renormalize momentum before the first step, and during each step.
            mom_norm = sqrt(yIn[3]*yIn[3] + yIn[4]*yIn[4] + yIn[5]*yIn[5]);
            for (int k = 3; k < 6; k ++){
               yIn[k] /= mom_norm;
            }
            
            myChawlaSharmaRKNstepperFullRenormalization->DumbStepper(yIn, step_len, yout);
            break;
         case 11:
            myNystromStepper->ComputeRightHandSide(yIn, dydx);
            myNystromStepper->Stepper(yIn,dydx,step_len,yout,yerr);
            break;
         default :
            myStepper->RightHandSide(yIn, dydx);                  
            //compute the value of dydx to supply to the stepper
            myStepper->Stepper(yIn,dydx,step_len,yout,yerr);   //call the stepper
      }
      
      //-> Then print the data
      // Modified to just print out X,Y,Z position data, 
      // so is easier to feed into Octave.
      for(int i=0; i<6;i++)
          if(columns[i])
              cout<<setw(20)<<yout[i]<<setw(20);
      
      //Copy yout into yIn
      for(int i=0; i<6;i++)
          yIn[i] = yout[i];
      cout<<"\n";
   }
   
   switch(stepper_no){
      case 7:
         delete myBogackiShampine23modified;
         break;
      case 8:
         delete myChawlaStepper;
         break;
      case 9:
         delete myChawlaStepper;
         break;
      case 10:
         delete myChawlaSharmaRKNstepperFullRenormalization;
         break;
      case 11:
         delete myNystromStepper;
         break;
      default:
         delete myStepper;
         break;
   }
    
/*-----------------END-STEPPING------------------*/    
}
