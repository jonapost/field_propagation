#include "G4UniformMagField.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4CashKarpRKF45.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "BogackiShampine23.hh"
#include "G4LineSection.hh"
#include "G4MagIntegratorStepper.hh"

//#include <system.h>
//#include "G4Types.h"

#define C 299.792458000

using namespace std;


//Made changes in yerr, heading not what you see in ouput


int main(int argc, char *args[]){

//	cout<<"#With Recommended Settings";

	int stepper_no = 2;
	G4double step_len = .1*CLHEP::mm;
	if(argc>1)
		stepper_no = atoi(args[1]);
	if(argc > 2)
		step_len = stof(args[2]);


	//1. Create a field :
	G4UniformMagField myField(0.01*tesla,0.0,0.0);

	//Create an Equation :
	G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);


    G4double particleCharge= +1.0*coulomb;  // in e+ units
    G4double spin=0.0;              // ignore the spin
    G4double magneticMoment= 0.0;   // ignore the magnetic moment

    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0);

	fEquation->SetChargeMomentumMass( chargeState,
	           0.01*kg*m*(1/s),
			1*gram);

	//Create a stepper :
	G4MagIntegratorStepper *myStepper, *exactStepper;


	//Choose the stepper based on the command line argument
	switch(stepper_no){
		case 1: myStepper = new G4CashKarpRKF45(fEquation);
			break;
		case 2: myStepper = new G4ExactHelixStepper(fEquation);
			break;
		case 3: myStepper = new BogackiShampine23(fEquation);
			break;
		default : myStepper = 0 ;
	}
    
    exactStepper = new G4ExactHelixStepper(fEquation);

	G4double yIn[] = {0.,0.,0.,0.01*kg*m*(1/s),0.,0.0};
    G4double dydx[] = {0,0,0,0,0,0};
	G4double yout[6], yerr[6];
    
//        -> First Print the (commented) title header
    
        cout<<"\n#  ";
        for(int i=0; i<6;i++){
            cout<<"yout["<<i<<"]"<<"\t\t";
        }

    
        //-> Then print the data
        cout<<"\n";
    
    cout.setf (ios_base::scientific);
    cout.precision(3);

        for(int j=0; j<25000; j++){
            
//            cout<<j+1<<" :\t\t";
            
        //compute the value of dydx to supply to the stepper :
        	myStepper->RightHandSide(yIn, dydx);
        //call the stepper
            myStepper->Stepper(yIn,dydx,step_len,yout,yerr);
            

        //-> Then print the data
        for(int i=0; i<6;i++){
            cout<<yout[i]<<"\t\t";
        }
        
        //Copy yout into yIn
            for(int i=0; i<6;i++){
                yIn[i] = yout[i];
            }
            
//
        cout<<"\n";
    }
    
    
    

    
    
    cout<<"\n\n#-------------End of output-----------------";

}
