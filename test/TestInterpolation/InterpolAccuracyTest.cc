
#include "G4UniformMagField.hh"         //For creating a uniform Magnetic Field
#include "G4QuadrupoleMagField.hh"
#include "G4CashKarpRKF45.hh"           //For ref. stepper
#include "G4ExactHelixStepper.hh"       //  - do -
#include "G4Mag_UsualEqRhs.hh"          //For equation object - fEquation

#include "G4MagIntegratorStepper.hh"    //For a stepper in general

#include "G4ExactHelixStepper.hh"       // Reference stepper

// #include "VernerRK56.hh"                //For Verner's RK 6(5) with interpolation
#include "DormandPrince745.hh"          //For Dormand-Prince's RK 5(4) Stepper with interpolation
#include "BogackiShampine45.hh"         //For Bogacki-Shampine's RK 5(4) Stepper with interpolation
// #include "G4BogackiShampine45.hh"         //For Bogacki-Shampine's RK 5(4) Stepper with interpolation
// #include "VernerRK78.hh"

#include <iomanip>

void     printout(G4double yout[], G4double youtX[], int columns[], G4double pass_no);
G4double toStepLen(G4double theta, G4double R);
G4double findRadius(G4double yOut2[],G4double yOut1[], G4double step_len);
void     test(G4double maxAngle, G4double Div, G4int verbose);

using namespace std;
using namespace CLHEP;

// Version 5.2 - Use Interpolation capability, and compare with exact solution

int main(int argc, char *argv[])
{
    G4int verbose = 1;
    
    cout<<"\t\t\t\t\t####### TEST FOR ACCURACY IN INTERPOLATION  ####### \n\n";
    
    for(int k=2; k< argc ;k++){
        test( atof(argv[1]), atof(argv[k]), verbose);
    }
    // Using the command line arguments for input :
    //    argv[1]              : maxAngle - The maximum rotation (e.g 2*pi for 1 full rotation)
    //    argv[2], argv[3] ... : Values of Div, step_len_in_theta = pi/Div

    //  test 
    //	  Last arg 0,1 : Verbosity, 0: result only, 1:Complete output of values
}

void printout(G4double yout[], G4double youtX[], int columns[], G4double pass_no)
{
    cout<<endl;
    //    cout.unsetf(ios_base::scientific);
    
    cout.setf (ios_base::scientific);
    cout.precision(4);
    cout<<setw(15)<<pass_no;
    
    for(int i=0; i<6;i++)
        if(columns[i]){
            cout<<setw(15)<<yout[i];
            cout<<setw(15)<<youtX[i];
            cout<<setw(15)<<yout[i] - youtX[i];
        }  
}

G4double toStepLen(G4double theta, G4double R = 0)
{
    if( R == 0.0 ){
       cerr<<"Error - R= 0.0";
    }
    
    return theta * R;
    
}

G4double findRadius(G4double yOut2[], G4double yOut1[], G4double step_len)
{
    G4double dy, dx, phi, R, theta;
    
    dx = yOut2[0] - yOut1[0];
    dy = yOut2[1] - yOut1[1];
    phi = atan2(dy, dx);
    theta = pi - 2.0*phi;
    R = step_len/theta;
    return R;
}

void test(G4double maxAngle, G4double Div, G4int verbose = 1)
{
    /* -----------------------------COMMON-------------------------------- */
    
    /****************************************************
     * This sections contains the initialisation of     *
     * various field and coordinate parameters :        *
     * initial position, momentum, Magnetic field,      *
     * and the declaration and definition of related    *
     * objects                                          *
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
    magneticMoment= 0.0,                // ignore the magnetic moment
    mass = 1.0;
    
    
    int
    columns[] =         //Variables in yOut[] you want to display - 0 for No, 1 for yes
    {
//   x_pos, y_pos, z_pos
//     1 , 1 , 1,                     
       0 , 0 , 0,
//   x_mom, y_mom, z_mom
//     0 , 0 , 0       
       1 , 1 , 1
    };
    
    //Initialising coordinates
    G4double yIn[7] = {x_pos, y_pos, z_pos, x_mom, y_mom, z_mom, 0.};
    G4double yInX[7] = {x_pos, y_pos, z_pos, x_mom, y_mom, z_mom, 0.};
    
    cout << "\n# Initial Vector yIn[] = ("
         << x_pos<<", "<<y_pos<<", "<<z_pos<<", "<<x_mom<<", "<<y_mom<<", "<<z_mom
         << ")" << endl;
    
    
    //Empty buckets for results
    G4double dydx[7] = {0., 0., 0., 0., 0., 0., 0.},
    dydxRef[7] = {0., 0., 0., 0., 0., 0., 0.},
    yout[7] = {0., 0., 0., 0., 0., 0., 0.},
    youtX[7] = {0., 0., 0., 0., 0., 0., 0.},
    yerr[7] = {0., 0., 0., 0., 0., 0., 0.},
    yerrX[7] = {0., 0., 0., 0., 0., 0., 0.},
    youtD[7] = {0., 0., 0., 0., 0., 0., 0.},
    youtDX[7] = {0., 0., 0., 0., 0., 0., 0.};
    
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

    
    /**-----------------------------  |     NEW     |  -------------------------------**/
    
    /****************************************************
     * Several calls to the stepper is made to find the *
     * values at the mesh points using a pre-defined    *
     * step length. Then more values are calculated at  *
     * the non-mesh points using interpolation.         *
     *                                                  *
     ****************************************************/
    
    /*-------PREP--------*/
    cout << "# Stepper information: " << endl;
    
    //Create a stepper :
    G4MagIntegratorStepper  *exactStepper;    // Reference stepper (for uniform field)
    
    // Creating the Reference solution ( Helix stepper )
    exactStepper = new G4ExactHelixStepper(fEquation);
    cout << "#      Reference solution: Exact Helix stepper " << endl;
    
    // The stepper to be tested
    // auto *myStepper = new G4BogackiShampine45(fEquation,6);
    // auto *myStepper = new BogackiShampine45(fEquation,6);
    // cout << "#      Trial solution: BogackiShampine45 Stepper" << endl;
    DormandPrince745 *myStepper = new DormandPrince745(fEquation, 6);
    cout << "#      Trial solution:     DormandPrince745 Stepper" << endl;
    //  auto *myStepper = new VernerRK78(fEquation,6);
    
    //  G4double *nextDydx = new G4double[7];
    G4double R;
    
    /* Determine the radius of curvature */
    G4double step_len = 20.0 *mm, yOut2[7];
    // myStepper->RightHandSide(yIn, dydx);
    // myStepper->Stepper(yIn,dydx,step_len,yOut2,yerr, nextDydx);
    
    exactStepper->RightHandSide(yIn, dydxRef);
    exactStepper->Stepper(yIn,dydxRef,step_len,yOut2,yerrX);
    R = findRadius(yOut2, yIn, step_len);
        
    /*---------xx--------*/

    /* NOW STEPPING */
    
    G4double theta = pi/Div;
    G4double angle = 0.;
    step_len = toStepLen(theta, R); //300.0 *mm;  //Step length in milimeters

    // G4double Err[6] = {0., 0., 0., 0., 0., 0.}; 
    G4double maxRelErr[6] =  {0., 0., 0., 0., 0., 0. };
    G4double maxAbsErr[6] =  {0., 0., 0., 0., 0., 0. };

    // Store values & parameters at maximum (absolute) difference
    G4double valueOutAbs[6]  =  {0., 0., 0., 0., 0., 0. };
    G4double valueRefAbs[6]  =  {0., 0., 0., 0., 0., 0. };
    G4double valAngleAbs[6]  =  {0., 0., 0., 0., 0., 0. };
    G4double valTauAbs[6]    =  {0., 0., 0., 0., 0., 0. };

    // Store values & parameters at maximum (relative) difference    
    G4double valueOutRel[6]  =  {0., 0., 0., 0., 0., 0. };
    G4double valueRefRel[6]  =  {0., 0., 0., 0., 0., 0. };
    G4double valAngleRel[6]  =  {0., 0., 0., 0., 0., 0. };
    G4double valTauRel[6]    =  {0., 0., 0., 0., 0., 0. };    
    
    G4double tau = 0.5, tau_step = 0.01*6.0/Div;      //0.01/Div
    
    //    G4double no_of_steps = 8 ;
    cout.unsetf(ios_base::scientific);

    cout<<"\n# Step_len = "<<step_len<<" mm "
    <<"\t theta_step = pi/"<<pi/theta;
    

    if(verbose>=1 ){
        //      - PRINTING HEADERS -
        cout<<"\n#";
        cout<<setw(13)<<"Rotation";
        for (int i=0; i<6;i++)
            if (columns[i])
                cout << setw(13)<< "yOut[" << i << "]"
                << setw(13)<< "yOutX["<< i << "]"
                << setw(13)<< "yOut-yOutX[" << i <<"]";
        
        cout<<"\n";
        cout<<"\n#  --------------- ";
        printout(yIn, yIn, columns, angle);
        cout<<"\n#  --------------- ";
    }

    
    myStepper->RightHandSide(yIn, dydx);
    exactStepper->RightHandSide(yInX, dydxRef);

    for(angle = theta; angle<maxAngle; angle += theta){
        
        //First calculate approximation at the mesh point
        
        myStepper->Stepper(yIn,dydx,step_len,yout,yerr);//, nextDydx);
        exactStepper->Stepper(yInX,dydxRef,step_len,youtX,yerrX);

        myStepper->SetupInterpolation();  // (yIn, dydx, step_len);
        //        Use a few interpolations in between the mesh points for dense output
        for(tau = tau_step; tau < 1.0; tau+=tau_step){

           myStepper->Interpolate( tau, youtD );
           //         ***********
            // Old Args: ( yIn, dydx, step_len, youtD, tau ); //Interpolating to tau
                                     
            exactStepper->Stepper(yInX, dydxRef, tau*step_len, youtDX, yerrX);
            
            //Checking if Error in this interpolated output is bigger
            for(int i =0; i<6; i++ )
            {
               // Relative error in Err[i]
               double relDiff= 0.0; // Err[i] = 0.0;            
               if(youtDX[i]!=0 && abs(youtD[i]/youtDX[i]) >= 1.0e-03){
                  //Neglecting values very close to zero - e.g. at multiples of pi
                  relDiff = (youtD[i] - youtDX[i])/youtDX[i];
               }
               // maxRelErr[i] = std::max( maxRelErr[i], std::fabs(relDiff) );
               if( std::fabs(relDiff) > maxRelErr[i] ){
                  maxRelErr[i] = std::fabs(relDiff);
                  // rawMaxErr[i] = relDiff;
                  valueOutRel[i]  = youtD[i];
                  valueRefRel[i]  = youtDX[i];
                  valAngleRel[i]  = angle;                  
                  valTauRel[i]    = tau;                  
               }
               
               double absDif= std::fabs(youtD[i] - youtDX[i] );
               // maxAbsErr[i] = std::max( maxAbsErr[i], absDif );                
               if( absDif > maxAbsErr[i] ) {
                  maxAbsErr[i] = absDif;
                  
                  valueOutAbs[i]  = youtD[i];
                  valueRefAbs[i]  = youtDX[i];
                  valAngleAbs[i]  = angle;                  
                  valTauAbs[i]    = tau;
               }               
            }
            if(verbose>=1)
            	printout(youtD, youtDX, columns, angle*tau);
        }
        
        if(verbose >=1){
            cout<<"\n#  --------------- ";
            printout(yout, youtX, columns, angle);
            cout<<"\n#  --------------- ";
        }
        
        for(int i=0; i<6; i++){
            yIn[i]  = yout[i];
            yInX[i] = youtX[i];
            // dydx[i] = nextDydx[i];
        }

        myStepper->RightHandSide(yIn, dydx);
        exactStepper->RightHandSide(yInX, dydxRef);
    }
    
    cout.unsetf(ios_base::scientific);
    
    cout<<"\n# \t\t\t  -------| RESULTS |--------\n";
    cout<<"\n# \tStep_len = "<<step_len<<" mm "
    <<"\t theta_step = pi/"<<pi/theta;
    
    cout.setf(ios_base::scientific);
    
    cout <<"\n# The Max. error in each coordinate : " << endl;
    cout << "#   Relative error: ";  
    for(int i=0;i <6; i++)
        cout << maxRelErr[i] << "  ";
    cout<<"\n";
    cout << "#  These occurred when : " << endl;    
    cout << "#   Reference value: ";      
    for(int i=0;i <6; i++)
        cout << valueOutRel[i] << " ";
    cout << endl;        
    cout << "#   New steper val : ";          
    for(int i=0;i <6; i++)
        cout << valueRefRel[i] << " ";
    cout << endl;

    // cout.unsetf(ios_base::scientific);
    cout << "#   Angle value    : ";              
    for(int i=0;i <6; i++)
        cout << valAngleRel[i] << " ";
    cout << endl;        
    cout << "#   Tau   value    : ";
    for(int i=0;i <6; i++)
        cout << valTauRel[i] << " ";
    cout << endl;  
    
    cout << endl;

    // cout.setf(ios_base::scientific);    
    cout << "#  Absolute error: ";  
    for(int i=0;i <6; i++)
        cout << maxAbsErr[i] << "  ";
    cout << endl;
    cout << "#  These occurred when : " << endl;
    cout << "#   Reference value: ";      
    for(int i=0;i <6; i++)
        cout << valueOutAbs[i] << " ";
    cout << endl;        
    cout << "#   New steper val : ";          
    for(int i=0;i <6; i++)
        cout << valueRefAbs[i] << " ";
    cout << endl;
    // cout.unsetf(ios_base::scientific);    
    cout << "#   Angle value    : ";              
    for(int i=0;i <6; i++)
        cout << valAngleAbs[i] << " ";
    cout << endl;        
    cout << "#   Tau   value    : ";
    for(int i=0;i <6; i++)
        cout << valTauAbs[i] << " ";
    cout << endl;            
    
    cout << "#\t\t\t\t\t\t------------- End of output ----------------- \n\n";
    
}
