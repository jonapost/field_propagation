//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//  Test for the accuracy of Interpolation method of RK Embedded steppers.
//
//  Usage:    InterpolAccuracyTest  max-angle  divisor-for-angle
//            InterpolAccuracyTest    2.5          4
//  Note:  This version includes check of result against
//
//    First version:  Somnath Banerjeee,  Summer 2015 -- GSoC 2015
//    Revisions:      John Apostolakis,  May 2016
//
//     - To compile use   "make G4TARGET=InterpolAccuracyTest"
//
//  Dev notes (open issues):     [ May 2016 ]
//     - Clarify role of 'max-angle' ie argv[0]
//     - Substitute calculation of radius (physics) for numerical estimation
// 
//  Future development:
//        Add argument parsing
//              -e value     --epsilon  value  --  relative error for step
//              -d value     --divisor  value  --  divisor for step 'angle'. angle = pi / divisor
//              -v           --verbose         --  verbose output
//              -N           --no-verbose      --  minimum output

#include "G4UniformMagField.hh"         //For creating a uniform Magnetic Field
#include "G4QuadrupoleMagField.hh"
#include "G4CashKarpRKF45.hh"           //For ref. stepper
#include "G4ExactHelixStepper.hh"       //  - do -
#include "G4Mag_UsualEqRhs.hh"          //For equation object - fEquation

#include "G4MagIntegratorStepper.hh"    //For a stepper in general

#include "G4ExactHelixStepper.hh"       // Reference stepper

#include "DormandPrince745.hh"          //For Dormand-Prince's RK 5(4) Stepper with interpolation
#include "G4BogackiShampine45.hh"         //For Bogacki-Shampine's RK 5(4) Stepper with interpolation

// #include "VernerRK56.hh"                //For Verner's RK 6(5) with interpolation
// #include "VernerRK78.hh"

#include <iomanip>

void     printout(G4double yout[], G4double youtX[], int columns[], G4double pass_no);
G4double toStepLen(G4double theta, G4double R);
G4double findRadius(G4double yOut2[],G4double yOut1[], G4double step_len);
void     test(G4double maxAngle, G4double Div, G4double epsilon, G4int verbose);

using namespace std;
using namespace CLHEP;

// Version 5.2 - Uses Interpolation capability, compare with exact solution   (S.B.)
// Version 5.3 - Print values for yOut (reference & trial stepper) at max difference (J.A.)
// Version 6.0 - Checks error in Step *before* testing interpolation          (J.A.)
//               ( Enables comparison of estimated step error with true interpolation error.)

int main(int argc, char *argv[])
{
    G4int verbose = 1;
    G4int startK = 2;
    double epsilon = 3.0e-4;  // It can be an arument (in future) e.g. = atof( argv[2] );

    /***
    if( argc > 3 ) { 
       double epsTry = atof( argv[2] );
       if( epsTry < 0.01 ) {
          epsilon = epsTry;
          cout << "# Using argv[2] = " << argv[2] << " as epsilon= " << epsTry << endl;
          startK ++;
       }
    } 
    
    if( argc > 4 ) {
       int verboseFlag= ( *argv[3] == 'v' || *argv[3] == 'V' );
       int noVerboseFlag= ( *argv[3] == 'n' || *argv[3] == 'N' );
       if( verboseFlag ) { verbose= 1; } 
       if( noVerboseFlag ) { verbose= 0; }
       if( verboseFlag || noVerboseFlag ) {
          startK++;
          cout << "# Step verbose = " << verbose << " from argv[3]= " << argv[3] << endl;
       }       
    }
     ***/
    
    cout<<"\t\t\t\t\t####### TEST FOR ACCURACY IN INTERPOLATION  ####### \n\n";
    
    for(int k= startK ; k< argc ;k++){
        test( atof(argv[1]), atof(argv[k]), epsilon, verbose);
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

void test(G4double maxAngle, G4double Div, G4double epsilon, G4int verbose = 1)
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
    
    // Arrays to store 'accumulants' - maximum, location/value, etc
    
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
    
    /****************************************************
     * Several calls to the stepper is made to find the *
     * values at the mesh points using a pre-defined    *
     * step length. Then more values are calculated at  *
     * the non-mesh points using interpolation.         *
     *                                                  *
     ****************************************************/
    
    /*------- Preparation --------*/
    cout << endl;
    cout << "# Stepper information: " << endl;
    cout << "#      Relative error 'epsilon' selected = " << epsilon << endl;

    //Create a stepper :
    G4MagIntegratorStepper  *exactStepper;    // Reference stepper (for uniform field)
    
    // Creating the Reference solution ( Helix stepper )
    exactStepper = new G4ExactHelixStepper(fEquation);
    cout << "#      Reference solution: Exact Helix stepper " << endl;
    
    // The stepper to be tested
    // auto *myStepper = new G4BogackiShampine45(fEquation,6);
    // cout << "#      Trial solution: BogackiShampine45 Stepper" << endl;
    DormandPrince745 *myStepper = new DormandPrince745(fEquation, 6);
    cout << "#      Trial solution:     DormandPrince745 Stepper" << endl;
    //  auto *myStepper = new VernerRK78(fEquation,6);
    
    //  G4double *nextDydx = new G4double[7];
    G4double R;
    
    /* Determine the radius of curvature */
    G4double step_len = 20.0 *mm;
    G4double yOut2[7];
    // myStepper->RightHandSide(yIn, dydx);
    // myStepper->Stepper(yIn,dydx,step_len,yOut2,yerr, nextDydx);
    
    exactStepper->RightHandSide(yIn, dydxRef);
    exactStepper->Stepper(yIn, dydxRef, step_len, yOut2, yerrX);
    R = findRadius(yOut2, yIn, step_len);
    cout << "# Radius found = " << R << endl;

    G4double tau = 0.5, tau_step = 0.01*6.0/Div;      //0.01/Div
    
    //    G4double no_of_steps = 8 ;
    cout.unsetf(ios_base::scientific);

    /* Parameters for steps */
    G4double angle = 0.;
    G4double theta = pi/Div;
    step_len = toStepLen(theta, R); //300.0 *mm;  //Step length in milimeters
    
    cout<<"\n# Step_len = " << step_len << " mm "
        <<"\t theta_step = pi/" << pi/theta << endl;

    if( verbose ){
        //    - Printing Header -
        cout << "\n#" << setw(13) << "Rotation";
        for (int i=0; i<6;i++)
            if (columns[i])
                cout << setw(13)<< "yOut[" << i << "]"
                     << setw(13)<< "yOutX["<< i << "]"
                     << setw(13)<< "yOut-yOutX[" << i <<"]";
        cout << endl << endl << "#  --------------- ";
        printout(yIn, yIn, columns, angle);
        cout << endl << "#  --------------- ";
    }

    myStepper->RightHandSide(yIn, dydx);
    exactStepper->RightHandSide(yInX, dydxRef);

    int stepNum = 0;
    for(angle = theta; angle<maxAngle; angle += theta)
    {
        stepNum++;

        // 1. Reference solution
        exactStepper->Stepper( yInX, dydxRef, step_len, youtX, yerrX );

        // 2. Calculate step - i.e. approximation at the mesh point
        myStepper->Stepper( yIn, dydx, step_len, yout, yerr ); //, nextDydx);

        // 3. Ensure that step succeeded - to given accuracy 'epsilon'
        double xErr= std::sqrt( yerr[0] * yerr[0] + yerr[1] * yerr[1] + yerr[2] * yerr[2] );
        double pErr= std::sqrt( yerr[3] * yerr[3] + yerr[4] * yerr[4] + yerr[5] * yerr[5] );
        double pMag= std::sqrt( yout[3] * yout[3] + yout[4] * yout[4] + yout[5] * yout[5] );

        bool  goodPos = xErr < epsilon * step_len;
        bool  goodMom = pErr < epsilon * pMag;
        bool goodStep= goodPos && goodMom;

        void ReportStepErrors( std::ostream &os, double xErr, double pErr, double pMag, double step_len);
        
        if( ! goodStep ) {
           cerr << "# -----------------------------------------------------------------------"
                << endl;
           cerr << "ERROR> Trial step error exceeds requested relative error limit = "
                << epsilon << endl;
           ReportStepErrors( std::cerr, xErr, pErr, pMag, step_len);
          /***           
           cerr << endl << endl;
                << endl;
           cerr << "         Error in Position:  Raw = " << xErr
                <<  " Relative = " << xErr / step_len << endl;
           cerr << "         Error in Momentum:  Raw = " << pErr
                <<  " Relative = " << pErr / pMag << endl;          
           cerr << "-----------------------------------------------------------------------"
                << endl; 
           ***/
           exit(1);
        }
        if( verbose || (stepNum==1) )
        {
           cout << "# -----------------------------------------------------------------------"
                << endl;
           cout << "#  Estimated integration error (obtained from stepper) in step " << stepNum << endl;
           ReportStepErrors( std::cout, xErr, pErr, pMag, step_len);
           /***
           cout << endl;
           cout << "#    Error in Position:  Raw = " << xErr
                <<  " Relative = " << xErr / step_len << endl;
           cout << "#    Error in Momentum:  Raw = " << pErr
                <<  " Relative = " << pErr / pMag << endl;            
           cout << "# -----------------------------------------------------------------------"
                << endl;
            ***/
        }        
        
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
    cout << "#  Raw error            : ";
    for(int i=0;i <6; i++)
        cout << maxAbsErr[i] << "  ";
    cout << endl;
    cout << "#  Scaled (relative) val: ";
    for(int i=0;i <6; i++) {
       // cout << maxAbsErr[i] / valueOutAbs[i] << "  ";       
       G4double ratio= 0.0;
       if( valueOutAbs[i] != 0.0 )
          ratio = maxAbsErr[i] / valueOutAbs[i];
       cout << ratio << "  ";
    }
    cout << endl;    
    cout << "#  These occurred when : " << endl;
    cout << "#   Reference value     : ";      
    for(int i=0;i <6; i++)
        cout << valueOutAbs[i] << " ";
    cout << endl;        
    cout << "#   New steper val      : ";          
    for(int i=0;i <6; i++)
        cout << valueRefAbs[i] << " ";
    cout << endl;
    // cout.unsetf(ios_base::scientific);    
    cout << "#   Angle value         : ";              
    for(int i=0;i <6; i++)
        cout << valAngleAbs[i] << " ";
    cout << endl;        
    cout << "#   Tau   value         : ";
    for(int i=0;i <6; i++)
        cout << valTauAbs[i] << " ";
    cout << endl;            
    
    cout << "#\t\t\t\t\t\t------------- End of output ----------------- \n\n";
    
}

void ReportStepErrors( std::ostream& os, double xErr, double pErr, double pMag, double step_len)
{
   // os << endl << endl;
   // os << "-----------------------------------------------------------------------"
   // << endl;
   os << "#         Error in Position:  Raw = " << xErr
      <<  " Relative = " << xErr / step_len << endl;
   os << "#         Error in Momentum:  Raw = " << pErr
      <<  " Relative = " << pErr / pMag << endl;          
   os << "# -----------------------------------------------------------------------"
      << endl;
}
