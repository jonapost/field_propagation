
#include "G4UniformMagField.hh"			//For creating a uniform Magnetic Field
#include "G4CashKarpRKF45.hh"			//For ref. stepper
#include "G4ExactHelixStepper.hh"		// 	- do -
#include "G4Mag_UsualEqRhs.hh"			//For equation object - fEquation
#include "VernerRK56.hh"				//For Verner's RK 6(5) with interpolation

#include "G4MagIntegratorStepper.hh"	//For a stepper in general

#include <iomanip>



using namespace std;
using namespace CLHEP;

//Version 4.0 - Test for Interpolation


int main(/*int argc, char *args[]*/){
    
    
    cout<<"\t\t\t\t\t\t\t####### TEST FOR INTERPOLATION  #######";
    
    /* -----------------------------COMMON-------------------------------- */
    
    /****************************************************
     * This sections contains the initialisation of 	*
     * various field and coordinate parameters :		*
     * initial position, momentum, Magnetic field,		*
     * and the declaration and definition of related 	*
     * objects											*
     ****************************************************/
    
    
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
    
    int
    columns[] =			//Variables in yOut[] you want to display - 0 for No, 1 for yes
    {
        1 ,  //x_pos
        1 ,  //y_pos
        1 ,  //z_pos
        0 ,  //x_mom
        0 ,  //y_mom
        0    //z_mom
        
        
    };
    

    
    //Initialising coordinates
    G4double yIn[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    G4double yInX[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    
    cout<<"\n# Initial Vector yIn[] = ("
        <<x_pos<<", "<<y_pos<<", "<<z_pos<<", "<<x_mom<<", "<<y_mom<<", "<<z_mom
        <<")";
    
    
    //Empty buckets for results
    G4double dydx[] = {0.,0.,0.,0.,0.,0.},
    dydxRef[] = {0.,0.,0.,0.,0.,0.},
    yout[6] = {0.,0.,0.,0.,0.,0.},
    youtX[6] = {0.,0.,0.,0.,0.,0.},
    yerr[6] = {0.,0.,0.,0.,0.,0.},
    yerrX[6] = {0.,0.,0.,0.,0.,0.};
    
    //1. Create a field :
    G4UniformMagField myField(G4ThreeVector(x_field, y_field, z_field));
    cout<<"\n# Field Info : "
    	<<"\n# \t Uniform Magnetic Field, (Bx, By, Bz) = "
    	<<"("<<x_field<<", "<< y_field<<", "<< z_field<<")";
    
    
    //Create an Equation :
    G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);
    
    
    
    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0);
    
    fEquation->SetChargeMomentumMass( chargeState,
                                     G4ThreeVector(x_mom, y_mom, z_mom).mag(), //momentum magnitude
                                     mass);//No place fo mass in fEquation though
    
    
    /***********************************************************************************************/
 
    

    
    
    
    
    
    
    /**-----------------------------  |     NEW		|  -------------------------------**/
    
    /****************************************************
     * This sections contains the creation and calls 	*
     * to the Stepper. Plus it also deals with printing	*
     * We do a single step with length specified by 	*
     * step_len. The approximations are done by calling *
     * the interpolate() function of the stepper.		*
     * The approximations at different fractions (tau)  *
     * of the step, as specified in the tau_list		*
     * The resulting values are compared with a ref.	*
     * stepper (like G4CashKarpRKF45) with repeated 	*
     * steps with reduces step sizes (step_len*tau)		*
     ****************************************************/
    
    
    
    			/*-------PREP--------*/
    G4double step_len = 10*mm;  //Step length in milimeters
    
    //Create a stepper :
    G4MagIntegratorStepper  *exactStepper;
//    G4MagIntegratorStepper *myStepper;
    
    //Creating the soulution stepper
    exactStepper = new G4CashKarpRKF45(fEquation);
    VernerRK56 *myStepper = new VernerRK56(fEquation);

    
                /*---------xx--------*/
    
    
    
    /* NOW STEPPING */

    G4double tau = 0.5;
    
    G4double tau_list[] = {0.0, 0.05, 0.113, 0.125, 0.339, 0.5, 0.775, 0.811, 0.9999, 1.0};
    
                /* PRINTING */
    
    //		- DETAILS -
    cout<<"\n# NOTE : The results with interpolation (in CERK Stepper) vs Reduced step sizes in Ref. Stepper";
    
    

    cout<<"\n#\n# Step Len = "<<step_len;
    
    //		- HEADERS -
    cout<<"\n#\n#";
    cout<<setw(13)<<"tau  ";
    for (int i=0; i<6;i++)
        if (columns[i])
            cout << setw(13)<< "yOut[" << i << "]"
            << setw(13) << "yOutX[" << i << "]"
            << setw(13) << "yOut-yOutX[" << i << "]";

    cout<<"\n";
    
    
    /* One step upto h=step_len */
    myStepper->RightHandSide(yIn, dydx);
    myStepper->Stepper(yIn,dydx,step_len,yout,yerr);

    /* ITERATING THROUGH tau_list */
    
    int j=0;
    
    while(j<10){
        

        tau = tau_list[j];
		
        /* Calling the ref. stepper */
        exactStepper->RightHandSide(yInX, dydxRef);
        exactStepper->Stepper(yInX,dydxRef,tau*step_len,youtX,yerrX);	//Stepping to a Reduced step size
        
        
        myStepper->interpolate(yIn, dydx, yout, step_len, tau);		//Interpolating to tau
    
        
        
   				 /*  - PRINTING -  */
        
        cout<<"\n";
        
        cout<<setw(13)<<tau;		//Send it below after setf, and precision calls to format it like others
        
        cout.setf (ios_base::scientific);
        cout.precision(3);
        for(int i=0; i<6;i++)
            if(columns[i]){
                cout<<setw(15)<<yout[i]<<setw(15);
                cout<<setw(15)<<youtX[i];
                cout<<setw(15)<<yout[i] - youtX[i];
            }
        cout.unsetf(ios_base::scientific);
        cout.precision(6);
        
        j++;

    }
    
    
    cout<<"\n\n#\t\t\t\t\t\t-------------End of output-----------------";
    
}
