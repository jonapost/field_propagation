
#include "G4UniformMagField.hh"			//For creating a uniform Magnetic Field
#include "G4QuadrupoleMagField.hh"
#include "G4CashKarpRKF45.hh"			//For ref. stepper
#include "G4ExactHelixStepper.hh"		// 	- do -
#include "G4Mag_UsualEqRhs.hh"			//For equation object - fEquation
#include "VernerRK56.hh"				//For Verner's RK 6(5) with interpolation
#include "DormandPrince745.hh" 			//For Dormand-Prince's RK 5(4) Stepper with interpolation
#include "BogackiShampine45.hh"			//For Bogacki-Shampine's RK 5(4) Stepper with interpolation
#include "G4MagIntegratorStepper.hh"	//For a stepper in general

#include <iomanip>


void printout(G4double yout[], int columns[], G4double pass_no);
G4double toStepLen(G4double theta, G4double R);
G4double findRadius(G4double yOut2[],G4double yOut1[], G4double step_len);


using namespace std;
using namespace CLHEP;

//Version 5.1 - Interpolation in Action (New : theta controlled)


int main(int argc, char *argv[]){
    
    
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
    
    x_field = 0.*gauss,               //Uniform Magnetic Field (x,y,z)
    y_field = 0.*gauss,
    z_field = 1000.*gauss ;
    
    
    //Set Charge etc.
    G4double particleCharge = +1.0,     // in e+ units
    spin=0.0,                           // ignore the spin
    magneticMoment= 0.0,                        // ignore the magnetic moment
    mass = 1.0;
    
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
    
    cout<<"\n# Initial Vector yIn[] = ("
    <<x_pos<<", "<<y_pos<<", "<<z_pos<<", "<<x_mom<<", "<<y_mom<<", "<<z_mom
    <<")";
    
    
    //Empty buckets for results
    G4double dydx[7] = {0.,0.,0.,0.,0.,0.};
    G4double yout[7] = {0.,0.,0.,0.,0.,0.};
    G4double youtD[7] = {0.,0.,0.,0.,0.,0.};
    
    G4double yerr[7] = {0.,0.,0.,0.,0.,0.};
    
    //1. Create a field :
    G4UniformMagField myField(G4ThreeVector(x_field, y_field, z_field));
    //    G4QuadrupoleMagField myField(z_field);
    
    
    cout<<"\n# Field Info : "
    <<"\n# \t Uniform Magnetic Field, (Bx, By, Bz) = "
    <<"("<<x_field<<", "<< y_field<<", "<< z_field<<")";
    //    cout<<"\n# Field Info : "
    //    <<"\n# \t Quadrupole Magnetic Field, Gradient = "
    //    << z_field;
    
    
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
     * Several calls to the stepper is made to find the *
     * values at the mesh points using a pre-defined 	*
     * step length. Then more values are calculated at	*
     * the non-mesh points using interpolation.			*
     * 													*
     ****************************************************/
    
    
    
    /*-------PREP--------*/
    
    //Create a stepper :
    //    G4MagIntegratorStepper  *exactStepper;
    //    G4MagIntegratorStepper *myStepper;
    
    //Creating the soulution stepper
    //    exactStepper = new G4CashKarpRKF45(fEquation);
    //    VernerRK56 *myStepper = new VernerRK56(fEquation,6);
    DormandPrince745 *myStepper = new DormandPrince745(fEquation, 6);
    G4double *nextDydx = new G4double[7];
    G4double R;
    
    /* Determine the radius of curvature */
    G4double step_len = 20.0 *mm, yOut2[7];
    myStepper->RightHandSide(yIn, dydx);
    myStepper->Stepper(yIn,dydx,step_len,yOut2,yerr, nextDydx);
    R = findRadius(yOut2, yIn, step_len);

    
    /*---------xx--------*/
    
    
    
    /* NOW STEPPING */
    G4double Div = 6.0;
    if(argc>1)
        Div = atoi(argv[1]);
    
    G4double theta = pi/Div;
    G4double angle = theta;
    step_len = toStepLen(theta, R); //300.0 *mm;  //Step length in milimeters
    
    G4double tau = 0.5;
    
//    G4double no_of_steps = 8 ;
    cout<<"\n# Step_len = "<<step_len<<" mm "
    	<<"\t theta_step = pi/"<<pi/theta;
    
    //		- PRINTING HEADERS -
    cout<<"\n#";
    cout<<setw(13)<<"Rotation";
    for (int i=0; i<6;i++)
        if (columns[i])
            cout << setw(13)<< "yOut[" << i << "]";
    
    cout<<"\n";
    
    myStepper->RightHandSide(yIn, dydx);
    for(angle = theta; angle<2.0*pi; angle += theta){
        
        //First calculate approximation at the mesh point
        
        myStepper->Stepper(yIn,dydx,step_len,yout,yerr, nextDydx);
        
        //Use a few interpolations in between the mesh points for dense output
        for(tau = 0.015; tau < 1.0; tau+=0.015){
            myStepper->interpolate(yIn,
                                   dydx,
                                   youtD,
                                   step_len,
                                   tau );			//Interpolating to tau
            printout(youtD, columns, angle*tau);
        }
        printout(yout, columns, angle);
        
        for(int i=0; i<6; i++){
            yIn[i] = yout[i];
            dydx[i] = nextDydx[i];
        }
        
    }
    
    cout<<"\n\n#\t\t\t\t\t\t-------------End of output-----------------";
    
}

void printout(G4double yout[], int columns[], G4double pass_no){
    
    cout<<endl;
//    cout.unsetf(ios_base::scientific);
    
    cout.setf (ios_base::scientific);
    cout.precision(4);
    cout<<setw(15)<<pass_no;

    
    for(int i=0; i<6;i++)
        if(columns[i]){
            cout<<setw(15)<<yout[i];
        }
    
}

G4double toStepLen(G4double theta, G4double R = 0){
    if(R == 0){
        cerr<<"Error";
        return 0;
    }
    
    return theta*R;
    
}

G4double findRadius(G4double yOut2[], G4double yOut1[], G4double step_len){
    G4double dy, dx, phi, R, theta;
    
    dx = yOut2[0] - yOut1[0];
    dy = yOut2[1] - yOut1[1];
    phi = atan(dy/dx);
    theta = pi - 2.0*phi;
    R = step_len/theta;
    return R;
}
