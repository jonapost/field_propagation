//
//  RKTest.cpp
//  Geant4
//
//  Created by Somnath on 14/07/15.
//
//

#include "RKTest.hh"
#include "RKTestDriver.cc"
//#include <iostream>
#include <iomanip>
using namespace CLHEP;

RKTest::RKTest() :
uField(G4UniformMagField(SetupUMF())),
QField(G4QuadrupoleMagField(SetupQMF())),
tTrack(G4FieldTrack(CreateTrack())){
//    fEqRhs = new G4Mag_UsualEqRhs(&uField);
    setEquation(&uField);
}

RKTest::~RKTest(){
    //Do nothing
}


void RKTest::testSteppersFixedQMF(string stepper_code, G4double step_len){
    tTrack = CreateTrack();
    QField = SetupQMF();
//    tDriver = SetDriver(stepper_code, QField, tTrack);
}


void RKTest::testSteppersFixed(string field_code = "umf",
                               string stepper_code = "ck45",
                               G4double step_len = 25.0*mm)
{
	//Code ..
}
void RKTest::print3(int columns[6],
                    int ifHeader,
                    G4double yOut[],
                    G4double yErr[],
                    G4double yOutX[]){
    
    if(ifHeader){
        cout<<"\n#";
        cout<<setw(7)<<"StepNo";
        for (int i=0; i<6;i++)
            if (columns[i])
                cout << setw(13)<< "yOut[" << i << "]"
                << setw(13) << "yErr[" << i << "]"
                << setw(13) << "yOut-yOutX[" << i << "]";
//        cout<<setw(13)<<"tan-1(y/x)";
        cout<<"\n";
    }
    
    else{
        cout.setf (ios_base::scientific);
        cout.precision(4);
        for(int i=0; i<6;i++)
            if(columns[i]){
                cout<<setw(15)<<yOut[i]<<setw(15);
                cout<<setw(15)<<yErr[i];
                cout<<setw(15)<<yOut[i] - yOutX[i];
            }
        cout.unsetf(ios_base::scientific);
        cout.precision(8);
        cout<<"\n";
    }
    
    
}

template < class STEPPER >
void RKTest::testSteppersFixedUMF( int columns[6],
                                   /*string stepper_code,*/
                                   G4double step_len ,
                                   int no_of_steps )
{
    
    tTrack = CreateTrack();
    setEquation(&uField);
    G4ExactHelixStepper *exactStepper = new G4ExactHelixStepper(fEqRhs);
    STEPPER *myStepper = new STEPPER(fEqRhs);
//    tDriver = SetDriver<STEPPER>(fEqRhs);//, tTrack);

    G4double
    yIn[12] = {0.,0.,0.,0.,0.,0.},
    yInX[12] = {0.,0.,0.,0.,0.,0.},
    dydx[12] = {0.,0.,0.,0.,0.,0.},
    dydxRef[12] = {0.,0.,0.,0.,0.,0.},
    yout[12] = {0.,0.,0.,0.,0.,0.},
    youtX[12] = {0.,0.,0.,0.,0.,0.},
    yerr[12] = {0.,0.,0.,0.,0.,0.},
    yerrX[12] = {0.,0.,0.,0.,0.,0.};
    

    tTrack.DumpToArray(yIn);
    tTrack.DumpToArray(yInX);

    //    -> First Print the (commented) title header
    
    print3(columns, 1, yout, yerr, youtX);

    
    /*----------------NOW STEPPING-----------------*/
    
    for(int j=0; j<no_of_steps; j++){
        
        cout<<setw(8)<<j + 1;           //Printing Step number
        
        myStepper->RightHandSide(yIn, dydx);    //compute the value of dydx to supply to the stepper
        
        myStepper->Stepper(yIn,dydx,step_len,yout,yerr);        //call the stepper
        
        
        exactStepper->RightHandSide(yInX, dydxRef);   //compute the value of dydx for the exact stepper
        
        exactStepper->Stepper(yInX,dydxRef,step_len,youtX,yerrX);       //call the exact stepper
        
        print3(columns, 0, yout, yerr, youtX);

        //Copy yout into yIn
        for(int i=0; i<12;i++){
            yIn[i] = yout[i];
            yInX[i] = youtX[i];
        }
        
        
    }
        
        
}

void RKTest::setEquation(G4MagneticField *pField){
    
    G4double momentum = tTrack.GetMomentum().mag();
    G4double
    particleCharge = 1.0,
    spin = 0.0,
    magneticMoment = 0.0 ;
    
    G4ChargeState chargestate(particleCharge, spin, magneticMoment);
    G4double mass = tTrack.GetRestMass();
    
//    if(fEqRhs)
//        delete fEqRhs;
//    fEqRhs = new G4Mag_UsualEqRhs(&uField);
    fEqRhs = new G4Mag_UsualEqRhs(pField);
    fEqRhs->SetChargeMomentumMass(chargestate, momentum, mass);
}

void RKTest::testPerformance(string stepper_code, string field_code){

    int max_no_loops = 20;
	G4double theEps = 1.;
    const G4double factor = 0.2;
    G4double physStep = 100.0*mm;
    G4int func_evals = 0;
        //    myDriver = SetupInt_Driver(stepper_no, fTrack );
        //testIntegrator(myDriver, fTrack, func_out[stepper_no]);
    cout<<"# ";
    cout<<setw(11)<<"-log10(eps)"
    	<<setw(15)<<"func_evals";
    cout<<"\n";
    for (int j=0; j<max_no_loops; j++) {

//        Reset();
        
        tTrack = CreateTrack();
        if(field_code == "qmf")
        	setEquation(&QField);
        else
            setEquation(&uField);
        

//        tDriver = SetDriver(stepper_code, fEqRhs, tTrack);
        tDriver->AccurateAdvance(tTrack, physStep, theEps,physStep); //Suggested h = physStep
        
        
        func_evals = tDriver->GetStepper()->GetfNoRHSCalls();

        cout<<setw(12)<< -log10(theEps);
        cout.setf (ios_base::scientific);
        cout.precision(4);
        cout<<setw(10)<<func_evals;
        cout.unsetf(ios_base::scientific);
        cout.precision(8);
        cout<<"\n";
        theEps *= factor;

    }
}

void RKTest::Reset(){
    uField = SetupUMF();
    QField = SetupQMF();

    if(MField)
        delete MField;
    tTrack = CreateTrack();
    if(tDriver)
        delete tDriver;
    if(fEqRhs)
        delete fEqRhs;
}

template <class STEPPER>
void RKTest::testAnyG4Stepper(string field_code){
    STEPPER *myStepper = new STEPPER(fEqRhs);
    G4int lalla = myStepper->GetfNoRHSCalls();
    cout<<"lalla = "<<lalla;
    
    int max_no_loops = 20;
    G4double theEps = 1.;
    const G4double factor = 0.2;
    G4double physStep = 100.0*mm;
    G4int func_evals = 0;
    //    myDriver = SetupInt_Driver(stepper_no, fTrack );
    //testIntegrator(myDriver, fTrack, func_out[stepper_no]);
    cout<<"# ";
    cout<<setw(11)<<"-log10(eps)"
    <<setw(15)<<"func_evals";
    cout<<"\n";
    for (int j=0; j<max_no_loops; j++) {
        
        //        Reset();
        
        tTrack = CreateTrack();
        if(field_code == "qmf")
            setEquation(&QField);
        else
            setEquation(&uField);
        
        
        tDriver = SetDriver<STEPPER>(fEqRhs);//, tTrack);
        tDriver->AccurateAdvance(tTrack, physStep, theEps,physStep); //Suggested h = physStep
        
        
        func_evals = tDriver->GetStepper()->GetfNoRHSCalls();
        
        cout<<setw(12)<< -log10(theEps);
        cout.setf (ios_base::scientific);
        cout.precision(4);
        cout<<setw(10)<<func_evals;
        cout.unsetf(ios_base::scientific);
        cout.precision(8);
        cout<<"\n";
        theEps *= factor;
    }
}

using namespace std;
int main(){
//    G4UniformMagField myField = SetupUMF();
//    G4FieldTrack myTrack = CreateTrack();
// 	G4MagInt_Driver *myDriver = SetDriver("ck45", myField, myTrack);
//    
//    cout<<myTrack.GetKineticEnergy();
    RKTest myTest;
    int columns[] = {1,0,0,0,0,0};
//    myTest.testSteppersFixedUMF(columns);
//    myTest.testPerformance("vr78");

//    G4CashKarpRKF45 *theStpr(0);
//    myTest.testAnyG4Stepper<G4CashKarpRKF45>("umf");
    myTest.testSteppersFixedUMF<G4CashKarpRKF45>(columns);
    
    cout<<"\n Hello";
    return 0;
}


