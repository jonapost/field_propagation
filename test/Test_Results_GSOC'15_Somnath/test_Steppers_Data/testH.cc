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
#include "DormandPrinceRK56.hh"

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
 5. G4ClassicalRK4
 6. G4SimpleHeum
 */


int main(int argc, char *args[]){
    
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* USER SPACE
     - Modify values here */
    
    int no_of_steps = 100;         //No. of Steps for the stepper
    int stepper_no = 7;             //Choose stepper no., for refernce see above
    G4double step_len = 50*mm;  //Step length in milimeters
    
    //Set coordinates here
    G4double
    x_pos = 0.,                 //pos - position
    y_pos = 0.,
    z_pos = 0.,
    
    x_mom = 0.,                 //mom - momentum
    y_mom = 10.,
    z_mom = 0.,
    
    x_field = 0.,               //Uniform Magnetic Field (x,y,z)
    y_field = 0.,
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
        step_len = (float)(stof(args[2])*mm);
    if(argc > 3)
        no_of_steps = atoi(args[3]);
    
    //Initialising coordinates
    G4double yIn[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    G4double yInX[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    
    
    //Empty buckets for results
    G4double dydx[7] = {0.,0.,0.,0.,0.,0.},
    dydxRef[7] = {0.,0.,0.,0.,0.,0.},
    yout[7] = {0.,0.,0.,0.,0.,0.},
    youtX[7] = {0.,0.,0.,0.,0.,0.},
    yerr[7] = {0.,0.,0.,0.,0.,0.},
    yerrX[7] = {0.,0.,0.,0.,0.,0.};
    
    //1. Create a field :
    G4UniformMagField myField(G4ThreeVector(x_field, y_field, z_field));
    
    //Create an Equation :
    G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);
    
    
    
    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0);
    
    fEquation->SetChargeMomentumMass( chargeState,
                                     G4ThreeVector(x_mom, y_mom, z_mom).mag(), //momentum magnitude
                                     mass);//No place fo mass in fEquation though
    
    //Create a stepper :
    G4MagIntegratorStepper *myStepper, *exactStepper;
    
    
    //Choose the stepper based on the command line argument
    switch(stepper_no){
        case 1: myStepper = new G4CashKarpRKF45(fEquation);
            break;
        case 0: myStepper = new G4ExactHelixStepper(fEquation);
            break;
//        case 2: myStepper = new BogackiShampine23(fEquation);
//            break;
//        case 3: myStepper = new DormandPrince745(fEquation);
//            break;
//        case 4: myStepper = new BogackiShampine45(fEquation);
//            break;
        case 5: myStepper = new G4ClassicalRK4(fEquation);
            break;
        case 6: myStepper = new G4SimpleHeum(fEquation);
            break;
        case 7: myStepper = new DormandPrinceRK56(fEquation);
            break;
        default : myStepper = 0 ;
    }
    
    //Creating the soulution stepper
    exactStepper = new G4ExactHelixStepper(fEquation);
    
    //        -> First Print the (commented) title header
    
    
    cout<<"\n#";
    cout<<setw(7)<<"StepNo";
    for (int i=0; i<6;i++)
        if (columns[i])
            cout << setw(13)<< "yOut[" << i << "]"
            << setw(13) << "yErr[" << i << "]"
            << setw(13) << "yOut-yOutX[" << i << "]";
    cout<<setw(13)<<"tan-1(y/x)";
    
    
    //-> Then print the data
    cout<<"\n";
    
    cout.setf (ios_base::scientific);
    cout.precision(3);
    
    /*-----------------------END PREPARING STEPPER---------------------------*/
    
    
    
    /*----------------NOW STEPPING-----------------*/
    
    for(int j=0; j<no_of_steps; j++){
        
        
        
        cout<<setw(8)<<j + 1;           //Printing Step number
        
        myStepper->RightHandSide(yIn, dydx);                    //compute the value of dydx to supply to the stepper
        myStepper->Stepper(yIn,dydx,step_len,yout,yerr);        //call the stepper
        
        
        exactStepper->RightHandSide(yInX, dydxRef);                 //compute the value of dydx for the exact stepper
        
        exactStepper->Stepper(yInX,dydxRef,step_len,youtX,yerrX);       //call the exact stepper
        
        
        
        //-> Then print the data
        cout.setf (ios_base::scientific);
        cout.precision(3);
        for(int i=0; i<6;i++)
            if(columns[i]){
                cout<<setw(15)<<yout[i]<<setw(15);
                cout<<setw(15)<<yerr[i];
                cout<<setw(15)<<yout[i] - youtX[i];
            }
        cout.unsetf(ios_base::scientific);
        cout.precision(6);
        cout<<setw(13)<<atan(yout[1]/yout[0])/degree;
        
        //Copy yout into yIn
        for(int i=0; i<6;i++){
            yIn[i] = yout[i];
            yInX[i] = youtX[i];
        }
        
        
        cout<<"\n";
    }
    
    /*-----------------END-STEPPING------------------*/
    
    
    
    
    
    
    cout<<"\n\n#-------------End of output-----------------";
    
}
