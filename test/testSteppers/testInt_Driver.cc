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



using namespace std;
using namespace CLHEP;

//Declare a bunch of functions
G4MagInt_Driver* SetupInt_Driver(G4int stepper_no,		//Create a new G4MagInt_Driver for us
                     			G4FieldTrack,
                     			G4UniformMagField);
void PrintResultOfStep(G4FieldTrack,		//Will print the (required) contents of G4FieldTrack
                       G4int stat_no);

void testIntegrator(G4MagInt_Driver*, G4FieldTrack&);			//The function for carrying out the test

G4FieldTrack CreateTrack();					//For creating a brand new G4FieldTrack

G4UniformMagField SetupField();


int main(int argc, char *argv[]){	//The main() function of the file

    
    //Some of the global variable could go here (?)
    
    G4int stepper_no = 2;		//Choosing a stepper_no, see above for reference

    G4FieldTrack fTrack = 0;
    G4UniformMagField myField = SetupField();
    
    G4MagInt_Driver *myDriver = SetupInt_Driver(stepper_no, fTrack, myField );
    fTrack = CreateTrack();
    
    if (argc>1) {					//Space for command line adjustments
        stepper_no = atoi(argv[1]);
    }
    
//    myDriver = SetupInt_Driver(stepper_no, fTrack );
    testIntegrator(myDriver, fTrack);
    
    //Printing of results - Done !
    return 0;
}




//Write the function definitions

void testIntegrator(G4MagInt_Driver *theDriver, G4FieldTrack &fTrack)
{
 
    
    G4int stat_no = 0;
    G4int no_of_steps = 10;
    G4double physStep = 50.*mm;
    G4double theEps = 1e-8;
    
    
    //Advance through the trajectory
    for(int i=0; i<no_of_steps; i++){
        theDriver->AccurateAdvance(fTrack, physStep, theEps,physStep/2.0); //Suggested h = physStep/2f
        PrintResultOfStep(fTrack, stat_no++);		//Print the resutls in each step
        
    }
    
    cout<<"\n# End of Output \n"
	    <<"\n# Actual no of steps done = "<<theDriver->GetNoTotalSteps();
    
    
    
    
    //G4Double hstep - the suggested step size
    //eps - allowed relative error
    //hinitial - initial suggested value of step, to be change automatically until the final pos reaches upto hstep advanced
    
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
    
    switch (Stepper_no) {
        case 1:
            pStepper = new G4ClassicalRK4(fEquation);
            break;
        case 2:
            pStepper = new DormandPrince745(fEquation);
            break;
        case 3:
            pStepper = new G4ExplicitEuler(fEquation);
            break;
        default:
            pStepper = 0;
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
    x_field = 0.,				//Uniform Magnetic Field (x,y,z)
    y_field = 0.,
    z_field = -0.1*tesla;
    
    
    
    G4UniformMagField myField(G4ThreeVector(x_field, y_field, z_field));
    return myField;
}

