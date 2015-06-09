
/* Some Unused headers present */

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
#include "G4MagIntegratorDriver.hh"
#include "G4PhysicalConstants.hh"
#include "G4ExplicitEuler.hh"



/*Add More Doc. */

/* What this is trying to do :
 
 Use various steppers, as them to give a certain accuracy
 Get the total no. of calls to the stepper
 Get the values for -
 A plot  of accuracy wanted vs total no of calls for stepper
 
 */

/* Stepper No.
 0. G4ExactHelix
 1. G4CashKarpRKF45
 2. BogackiShampine23
 3. DormandPrince745
 4. BogackiShampine45
 5. G4ClassicalRK4
 6. G4SimpleHeum
 7. G4ExplicitEuler
 */



#define G4ExactHelix_stages 1
#define G4CashKarpRKF45_stages 6
#define BogackiShampine23_stages 4
#define DormandPrince745_stages 7
#define BogackiShampine45_stages 8
#define G4ClassicalRK4_stages 4
#define G4SimpleHeum_stages 2
#define G4ExplicitEuler_stages 1

G4int stages;
G4String stepper_name;

using namespace std;
using namespace CLHEP;

//Declare a bunch of functions
G4MagInt_Driver* SetupInt_Driver(G4int stepper_no,		//Create a new G4MagInt_Driver for us
                                 G4FieldTrack,
                                 G4UniformMagField);
void PrintResultOfStep(G4FieldTrack,		//Will print the (required) contents of G4FieldTrack
                       G4int stat_no);

void testIntegrator(/*G4MagInt_Driver*, G4FieldTrack&, G4int* */);			//The function for carrying out the test

G4FieldTrack CreateTrack();					//For creating a brand new G4FieldTrack

G4UniformMagField SetupField();


int main(/*int argc, char *argv[]*/){	//The main() function of the file
    
    testIntegrator();
    
    return 0;
}




//Write the function definitions

void testIntegrator(/*G4MagInt_Driver *theDriver, G4FieldTrack &fTrack, G4int fout[]*/)
{
    
    G4int func_evals = 0;
    G4double physStep = 2000.*mm;
    G4double theEps = 1.;
    
    
    G4MagInt_Driver *theDriver = 0;
    G4UniformMagField myField = SetupField();
    G4FieldTrack fTrack = 0;
    G4int fout[8][20];//   [8][10];		//fout [max(stepper_no)] [no_of_loops]
    G4double epsOut[20];
    
    const G4int
    max_no_steppers = 7,
    max_no_loops = 20;
    const G4double factor = 0.2;
    //    fout = new G4int[max_no_steppers][max_no_loops];
    
    
    
    cout<<"# Test for no. of function evaluation for a given accuracy";
    cout<<"\n# Curve length = "<<physStep<<endl;
    
    cout<<"\n#";
    cout<<setw(10)<<"-log(eps)";
    for(int stepper_no=0; stepper_no<=max_no_steppers; stepper_no++)
    {
        theEps = 1.;
        
        //    myDriver = SetupInt_Driver(stepper_no, fTrack );
        //testIntegrator(myDriver, fTrack, func_out[stepper_no]);
        for (int j=0; j<max_no_loops; j++) {
            
            
            epsOut[j] = -log10(theEps);	//Evaluated too many times, assigned same set of values for all stepper_no : 0:7
            fTrack = CreateTrack();
            theDriver = SetupInt_Driver(stepper_no, fTrack, myField );
            theDriver->AccurateAdvance(fTrack, physStep, theEps,physStep); //Suggested h = physStep
            
            func_evals = theDriver->GetTotalNoStepperCalls();//stages*theDriver->GetTotalNoStepperCalls() + theDriver->GetNoTotalSteps();
            fout[stepper_no][j] = func_evals;
            theEps *= factor;
            
            
        }
        cout<<setw(19)<<stepper_name;
    }
    
    cout<<endl;
    
    for(int j=0; j<max_no_loops; j++){
        cout<<setw(10)<<epsOut[j];
        for(int k=0; k<=max_no_steppers; k++){
            cout<<setw(19)<<fout[k][j];
        }
        cout<<endl;
    }
    
}


//Needs improvement though
G4MagInt_Driver* SetupInt_Driver(G4int Stepper_no, G4FieldTrack theTrack, G4UniformMagField myField){
    //Create a field :
    
    //= new G4UniformMagField (G4ThreeVector(x_field, y_field, z_field));
    
    G4double particleCharge = +1.0, 	// in e+ units
    spin=0.0,              						// ignore the spin
    magneticMoment= 0.0,   						// ignore the magnetic moment
    mass = theTrack.GetRestMass();//proton_mass_c2;					//1.67262178e-27*kilogram*c*c;
    //create a Chargestate
    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin,
                              magneticMoment);
    
    
    G4double hminimum = 0.001*mm;
    G4int
    numComponents = 6,
    statisticsVerbose = 0;
    G4double momentum = theTrack.GetMomentum().mag();// Needs reverting
    
    //Create an Equation :
    G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);
    
    
    
    fEquation->SetChargeMomentumMass ( chargeState,
                                      momentum, //momentum magnitude
                                      mass);//No place fo mass in fEquation though
    
    
    G4MagIntegratorStepper *pStepper;
    
    switch(Stepper_no){
        case 1: pStepper = new G4CashKarpRKF45(fEquation);
            stages = G4CashKarpRKF45_stages-1;
            stepper_name = "G4CashKarpRKF45";
            break;
        case 0: pStepper = new G4ExactHelixStepper(fEquation);
            stages = G4ExactHelix_stages - 1;
            stepper_name = "Exact";
            break;
        case 2: pStepper = new BogackiShampine23(fEquation);
            stages = BogackiShampine23_stages - 1;
            stepper_name = "BogackiShampine23";
            break;
        case 3: pStepper = new DormandPrince745(fEquation);
            stages = DormandPrince745_stages - 1;
            stepper_name = "DormandPrince745";
            break;
        case 4: pStepper = new BogackiShampine45(fEquation);
            stages = BogackiShampine45_stages - 1;
            stepper_name = "BogackiShampine45";
            break;
        case 5: pStepper = new G4ClassicalRK4(fEquation);
            stages = G4ClassicalRK4_stages - 1 ;
            stepper_name = "G4ClassicalRK4";
            break;
        case 6: pStepper = new G4SimpleHeum(fEquation);
            stages = G4SimpleHeum_stages - 1;
            stepper_name = "G4SimpleHeum";
            break;
        case 7: pStepper = new G4ExplicitEuler(fEquation);
            stages = G4ExplicitEuler_stages - 1;
            stepper_name = "G4ExplicitEuler";
            break;
        default : pStepper = 0 ;
            G4cerr<<" Invalid Stepper_no";
            break;
    }
    
    //Setup the Driver
    G4MagInt_Driver *theDriver = new G4MagInt_Driver(hminimum,
                                                     pStepper,
                                                     numComponents,
                                                     statisticsVerbose);
    
    
    return theDriver;
}

void PrintResultOfStep( G4FieldTrack theTrack, G4int stat_no/* Pass in the fTrack*/){
    //Call the Printstatus method of the Driver class (?)
    
    //Have an option for file output (?)
    
    G4ThreeVector
    position = theTrack.GetPosition(),
    momemntum = theTrack.GetMomentum();
    
    
    
    G4double yOut[24];
    
    //    static int stat_no = 0;
    
    if(stat_no++ == 0){
        //Print the headers
        cout
        << std::setw( 12) << "X(mm)" << " "
        << std::setw( 12) << "Y(mm)" << " "
        << std::setw( 12) << "Z(mm)" << " "
        << std::setw( 12) << " P_x " << " "
        << std::setw( 12) << " P_y " << " "
        << std::setw( 12) << " P_z " << " ";
        
        
        cout<<endl;
        
    }
    
    cout.setf (ios_base::scientific);
    cout.precision(3);
    
    theTrack.DumpToArray(yOut);
    
    cout<<setw(13)<<position.x();
    cout<<setw(13)<<position.y();
    cout<<setw(13)<<position.z();
    cout<<setw(13)<<momemntum.x();
    cout<<setw(13)<<momemntum.y();
    cout<<setw(13)<<momemntum.z();
    
    
    cout<<endl;
    
}
G4FieldTrack CreateTrack(){
    
    
    //Setup a field
    G4double
    x_mom = 0.,					//mom - momentum
    y_mom = 10.,//0.5*proton_mass_c2,
    z_mom = 0.;
    
    G4double mass = 1;
    
    
    G4double momentum = G4ThreeVector(x_mom,y_mom, z_mom).mag();
    //for the G4FieldTrack - momentum, position etc.
    
    G4ThreeVector position(0.,0.,0.);
    
    G4ThreeVector UnitMomentum = G4ThreeVector(x_mom,y_mom, z_mom).unit();
    
    //    G4double mass = 1.67262178e-27*kilogram;
    G4double kineticEnergy = sqrt(sqr(momentum) + sqr(mass)) - mass;
    G4double velocity = momentum/(mass + kineticEnergy);
    G4double labTof = 10.0*ns ;
    G4double properTof = 0.1*ns;
    
    
    
    
    //Declare and construct the fieldtrack :
    G4FieldTrack f( position,
                   UnitMomentum,
                   0.0,	//for S curve len - starting point
                   kineticEnergy,
                   mass,
                   velocity,
                   labTof,
                   properTof,
                   0	//for Spin
                   );
    
    return f;
    
    
}


G4UniformMagField SetupField(){
    
    G4double
    x_field = 0.1*tesla,				//Uniform Magnetic Field (x,y,z)
    y_field = 0.1*tesla,
    z_field = -0.1*tesla;
    
    
    
    G4UniformMagField myField(G4ThreeVector(x_field, y_field, z_field));
    return myField;
}








