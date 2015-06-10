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
#include "ChawlaSharmaRKNstepper_with_renormalization.hh"
#include "modifiedBogackiShampine23.hh"
#include "G4NystromRK4.hh"

//#include
//#include <system.h>
//#include "G4Types.h"

using namespace std;
using namespace CLHEP;

//Version 2.0 - Includes direct comparison with G4ExactHelix

/* Stepper No.
    0. G4ExactHelix
    1. G4CashKarpRKF45
    2. BogackiShampine23
    3. DormandPrince745
    4. BogackiShampine45
 */


int main(int argc, char *args[]){

    /* -----------------------------SETTINGS-------------------------------- */
    
/* USER SPACE
    - Modify values here */
    
    int no_of_steps = 1000;         //No. of Steps for the stepper
    int stepper_no = 3;             //Choose stepper no., for refernce see above
    G4double step_len = 10*mm;  //Step length in milimeters
  
    //Set coordinates here
    G4double
    x_pos = 0.,                 //pos - position
    y_pos = 0.,
    z_pos = 0.,

    x_mom = 4.,                 //mom - momentum
    y_mom = 10.,
    z_mom = 2.,
    
    x_field = 0.05*tesla,               //Uniform Magnetic Field (x,y,z)
    y_field = 0.2*tesla,
    z_field = -0.1*tesla ;
    
    
    //Set Charge etc.
    G4double particleCharge = +1.0,     // in e+ units
    spin=0.0,                                   // ignore the spin
    magneticMoment= 0.0,                        // ignore the magnetic moment
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
    cout << "Test -4 \n";
    //Initialising coordinates
    G4double yIn[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    G4double yInX[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    
    
    //Empty buckets for results
    G4double dydx[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        dydxRef[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        yout[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        youtX[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        yerr[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        yerrX[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

    //1. Create a field :
    G4UniformMagField myField(G4ThreeVector(x_field, y_field, z_field));

    //Create an Equation :
    G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);



    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0);

    fEquation->SetChargeMomentumMass( chargeState,
                                G4ThreeVector(x_field, y_field, z_field).mag(), //momentum magnitude
                                mass);//No place fo mass in fEquation though
    
    //Create a stepper :
    G4MagIntegratorStepper *myStepper, *exactStepper;
    ChawlaSharmaRKNstepper *myChawlaStepper;
    ChawlaSharmaRKNstepper_with_renormalization *myChawlaStepper_with_renormalization;
    modifiedBogackiShampine23 *myModifiedBogackiShampine23;
    cout << "Test -3 \n";
    G4NystromRK4 *myNystromStepper;
    cout << "Test -2 \n";

    //Choose the stepper based on the command line argument
    switch(stepper_no){
        case 1: myStepper = new G4CashKarpRKF45(fEquation);
            break;
        case 0: myStepper = new G4ExactHelixStepper(fEquation);
            break;
        case 2: myStepper = new BogackiShampine23(fEquation);
            break;
        case 3: myStepper = new DormandPrince745(fEquation);
            break;
        case 4: myStepper = new BogackiShampine45(fEquation);
            break;
        case 5: myStepper = new G4ClassicalRK4(fEquation);
            break;
        case 6: myStepper = new G4SimpleHeum(fEquation);
            break;
        case 7: myModifiedBogackiShampine23 = new modifiedBogackiShampine23(fEquation); // No renormalization inside of a step.
            break;
        case 8: myChawlaStepper = new ChawlaSharmaRKNstepper(fEquation); // No renormalization of momentum.
            break;
        case 9: myChawlaStepper = new ChawlaSharmaRKNstepper(fEquation); // Renormalization in between steps.
            break;
        case 10: 
                myChawlaStepper_with_renormalization = new ChawlaSharmaRKNstepper_with_renormalization(fEquation);
                // Normalization inbetween steps, and inside of steps.
            break;
        case 11:
                cout << "Test -1 \n";
                myNystromStepper = new G4NystromRK4(fEquation);
                cout << "Test 0 \n";
            break;
        default : myStepper = 0 ;
    }
    
    //Creating the soulution stepper
    exactStepper = new G4ExactHelixStepper(fEquation);
    
//        -> First Print the (commented) title header
    
    /*
        cout<<"\n#";
        cout<<setw(7)<<"StepNo";
    for (int i=0; i<6;i++)
        if (columns[i])
            cout << setw(13)<< "yOut[" << i << "]"
                 << setw(13) << "yErr[" << i << "]"
                 << setw(13) << "yOut-yOutX[" << i << "]";

    
        //-> Then print the data
        cout<<"\n";
    */
    //cout.setf (ios_base::scientific);
    //cout.precision(5);
    
    // For output to Octave (for visualization)
    cout.setf (ios_base::fixed);
    cout.precision(9);
    
/*-----------------------END PREPARING STEPPER---------------------------*/
    

    
    /*----------------NOW STEPPING-----------------*/
    
    cout << " Test 1 \n";
    
        G4double mom_norm;
        
        for(int j=0; j<no_of_steps; j++){
        
            //cout<<setw(8)<<j + 1;           //Printing Step number
            switch(stepper_no){
                case 7:
                  //myModifiedBogackiShampine23->modifiedRightHandSide(yIn, dydx); // Temp change to test a theory
                  myModifiedBogackiShampine23->RightHandSide(yIn, dydx); // Renormalize between steps ( Standard RightHandSide() )
                  myModifiedBogackiShampine23->Stepper(yIn, dydx, step_len, yout, yerr);
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
                  
                  myChawlaStepper_with_renormalization->DumbStepper(yIn, step_len, yout);
                  break;
               
               case 11:
                  cout << "Test 2 \n";
                  myNystromStepper->ComputeRightHandSide(yIn, dydx);
                  cout << "Test 3 \n";
                  myNystromStepper->Stepper(yIn,dydx,step_len,yout,yerr);
                  cout << "Test 4 \n";
                  break;
                  
               default :
                  myStepper->RightHandSide(yIn, dydx);                    //compute the value of dydx to supply to the stepper
                  myStepper->Stepper(yIn,dydx,step_len,yout,yerr);        //call the stepper
            }
            
            // For these tests I don't really care about the exact stepper
            // because I'm comparing it to BogackiShampine23
            //exactStepper->RightHandSide(yInX, dydxRef);                 //compute the value of dydx for the exact stepper
                                                                            // ^ Same dydx space being used
            //exactStepper->Stepper(yInX,dydxRef,step_len,youtX,yerrX);       //call the exact stepper
            
            
            
        //-> Then print the data
        // Modified to just print out X,Y,Z position data, so is easier to feed into Octave.
            for(int i=0; i<6;i++)
                if(columns[i]){
                    cout<<setw(20)<<yout[i]<<setw(20);
                    //cout<<setw(15)<<yerr[i];
                    //cout<<setw(15)<<yout[i] - youtX[i];
                }
        
        //Copy yout into yIn
            for(int i=0; i<6;i++){
                yIn[i] = yout[i];
                yInX[i] = youtX[i];
            }
            

        cout<<"\n";
        }
    
    /*-----------------END-STEPPING------------------*/
    
    
    

    
    
    //cout<<"\n\n#-------------End of output-----------------";

}
