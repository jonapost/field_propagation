//
//  RKTest.h
//  Geant4
//
//  Created by Somnath on 14/07/15.
//

#ifndef __Geant4__RKTest__
#define __Geant4__RKTest__
/*
 Here we need different types of tests to be run : 
 1. Plain test of the stepper : Takes the name of the stepper and gives
 	the 3 line output per positive choice of a coordinate
*/

#include "RKTestField.hh"
#include "RKTestTrack.hh"
#include "RKTestDriver.hh"
#include "../src/RKTestDriver.cc"
#include <stdio.h>
#include <iomanip>

using namespace std;
using namespace CLHEP;

class RKTest{
private:
    //Members - all needed for the program Run
//    G4UniformMagField uField ;
//    G4QuadrupoleMagField QField ;
    G4MagneticField *MField;
    G4FieldTrack tTrack ;
    G4MagInt_Driver *tDriver ;
    FSALMagInt_Driver *FDriver;
    G4Mag_UsualEqRhs *fEqRhs;
public:
    //Member functions
    RKTest();
    ~RKTest();
    int RunTest();	//Run the test and return a value indicator
    void Reset();
    void setEquation(G4MagneticField *pField);
    

//    
    static void print3(int columns[], int ifHeader , G4double yOut[], G4double yErr[], G4double yOutX[]);
//    void print3Data(G4double yOut[], yErr[], yOutX);
    
    template < class STEPPER >
    void testSteppersFixedUMF(int columns[6],
                              /*string stepper_code = "ck45",*/
                              G4double factor = 1.,
                              G4double step_len = 25.0*CLHEP::mm,
                              int no_of_steps = 100)
    {
        
        tTrack = CreateTrack();
        MField = new G4UniformMagField(SetupUMF(factor));
        setEquation(MField);
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
    
    template<class STEPPER>
    void testG4Stepper(std::string field_code = "umf", G4double factor = 1.)
    {
        
        //    STEPPER *myStepper = new STEPPER(fEqRhs);
        
        int max_no_loops = 20;
        G4double theEps = 1.;
        const G4double accu_factor = 0.2;
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
            if(field_code == "qmf"){
//                QField = SetupQMF(factor);
                MField = new G4QuadrupoleMagField(SetupQMF(factor));
//                setEquation(&QField);
                setEquation(MField);
                
            }
            else{
//                uField = SetupUMF(factor);
                MField = new G4UniformMagField(SetupUMF(factor));
//                setEquation(&uField);
                setEquation(MField);
            }

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
            theEps *= accu_factor;
        }
    }
    
    template<class STEPPER>
    void testFSALStepper(std::string field_code = "umf", G4double factor = 1.)
    {
        
        //    STEPPER *myStepper = new STEPPER(fEqRhs);
        
        int max_no_loops = 20;
        G4double theEps = 1.;
        G4double accu_factor = 0.2;
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
            if(field_code == "qmf"){
                //                QField = SetupQMF(factor);
                MField = new G4QuadrupoleMagField(SetupQMF(factor));
                //                setEquation(&QField);
                setEquation(MField);
                
            }
            else{
                //                uField = SetupUMF(factor);
                MField = new G4UniformMagField(SetupUMF(factor));
                //                setEquation(&uField);
                setEquation(MField);
            }
            
            /*****************	Only Difference from testG4Driver ********************/
            FDriver = SetFSALDriver<STEPPER>(fEqRhs);//, tTrack);
            FDriver->AccurateAdvance(tTrack, physStep, theEps,physStep); //Suggested h = physStep
            /*************************************************************************/
            
            func_evals = FDriver->GetStepper()->GetfNoRHSCalls();
            
            cout<<setw(12)<< -log10(theEps);
            cout.setf (ios_base::scientific);
            cout.precision(4);
            cout<<setw(10)<<func_evals;
            cout.unsetf(ios_base::scientific);
            cout.precision(8);
            cout<<"\n";
            theEps *= accu_factor;
        }
    }

    template<class STEPPER, class REF_STEPPER>
    void testStepperInterpolant(int columns[6],
                                G4double step_len_pi_divisor = 6.0,
                                G4double maxAngle = 2.0*CLHEP::pi,
                                G4double factor = 1.,
                                std::string field_code = "umf")
    {
        bool ifPrintany = false;
        
        tTrack = CreateTrack();
        if(field_code == "qmf"){
            //                QField = SetupQMF(factor);
            MField = new G4QuadrupoleMagField(SetupQMF(factor));
            //                setEquation(&QField);
            setEquation(MField);
            
        }
        else{
            //                uField = SetupUMF(factor);
            MField = new G4UniformMagField(SetupUMF(factor));
            //                setEquation(&uField);
            setEquation(MField);
        }
    
        
        REF_STEPPER *exactStepper = new REF_STEPPER(fEqRhs);
        STEPPER *myStepper = new STEPPER(fEqRhs);
        
        G4double
        yIn[12] = {0.,0.,0.,0.,0.,0.},
        yInX[12] = {0.,0.,0.,0.,0.,0.},
        dydx[12] = {0.,0.,0.,0.,0.,0.},
        dydxRef[12] = {0.,0.,0.,0.,0.,0.},
        yout[12] = {0.,0.,0.,0.,0.,0.},
        youtX[12] = {0.,0.,0.,0.,0.,0.},
        yerr[12] = {0.,0.,0.,0.,0.,0.},
        yerrX[12] = {0.,0.,0.,0.,0.,0.},
        youtD[7] = {0., 0., 0., 0., 0., 0., 0.},
        youtDX[7] = {0., 0., 0., 0., 0., 0., 0.};
        
        
        tTrack.DumpToArray(yIn);
        tTrack.DumpToArray(yInX);
        
        
        /* Determine the radius of curvature */
        G4double step_len = 20.0 *mm, yOut2[7];
        // myStepper->RightHandSide(yIn, dydx);
        // myStepper->Stepper(yIn,dydx,step_len,yOut2,yerr, nextDydx);
        
        exactStepper->RightHandSide(yIn, dydxRef);
        exactStepper->Stepper(yIn,dydxRef,step_len,yOut2,yerrX);
        G4double R = findRFrom2pt(yOut2, yIn, step_len);
        
        G4double theta = pi/step_len_pi_divisor;
        G4double angle = 0.;
        step_len = toStepLen(theta, R); //300.0 *mm;  //Step length in milimeters
        
        //Printing the headers
        for(int i=0; i<6; i++){
            if(columns[i]){
                ifPrintany = true;
                cout<<"# "
                <<setw(10)<<"yOut["<<i<<"]"
                <<setw(14)<<"yOutX["<<i<<"]"
                <<setw(16)<<"yOut - yOutX["<<i<<"]";
            }
        }
        
        myStepper->RightHandSide(yIn, dydx);
        exactStepper->RightHandSide(yInX, dydxRef);

        if(ifPrintany)
            print3(columns, 0, yout, yerr, yout);
        
        
        G4double Err[6] = {0., 0., 0., 0., 0., 0.},
        maxErr[6] =  {0., 0., 0., 0., 0., 0. },
        avgErr[6] = {0., 0., 0., 0., 0., 0. };
        G4double tau = 0.5, tau_step = 0.01*6.0/step_len_pi_divisor;//0.01/Div
        int no_of_steps = 0;
        
        for(angle = theta; angle<maxAngle; angle += theta){
            
            //First calculate approximation at the mesh point
            
            myStepper->Stepper(yIn,dydx,step_len,yout,yerr);//, nextDydx);
            exactStepper->Stepper(yInX,dydxRef,step_len,youtX,yerrX);
//            myStepper->SetupInterpolate_high(yIn, dydx, step_len);
            myStepper->SetupInterpolate(yIn, dydx, step_len);	//--Original
            //        Use a few interpolations in between the mesh points for dense output
            for(tau = tau_step; tau < 1.0; tau+=tau_step){
                no_of_steps++;
//                myStepper->Interpolate_high(yIn,
//                                       dydx,
//                                       step_len,
//                                       youtD,
//                                       tau );           //Interpolating to tau
                myStepper->Interpolate(yIn, //--Original
                                            dydx,
                                            step_len,
                                            youtD,
                                            tau );           //Interpolating to tau
                
                exactStepper->Stepper(yInX,dydxRef,tau*step_len,youtDX,yerrX);
                
                //Checking if Error in this interpolated output is bigger
                for(int i =0; i<6; i++ )
                {
                    if(youtDX[i]==0 || abs(youtDX[i]/youtD[i])< 1.0e-03){
                        //Neglecting values ridiculously close to zero at multiples of pi
                        Err[i] = 0.0;
                        continue;
                        
                    }
                    Err[i] = (youtD[i] - youtDX[i])/youtDX[i];
                    Err[i] = abs(Err[i]);
                    if(Err[i] > maxErr[i]) maxErr[i] = Err[i] ;
                    avgErr[i] += Err[i];
                    
                }
                if(ifPrintany)
                	print3(columns, 0, youtD, youtDX, youtDX);//(youtD, youtDX, columns, angle*tau);
            }
            if(ifPrintany){
                cout<<"\n#  --------------- \n";
                print3(columns, 0, yout, youtX, youtX);
                cout<<"\n#  --------------- \n";
            }

            
            for(int i=0; i<6; i++){
                yIn[i] = yout[i];
                yInX[i] = youtX[i];
            }
            
            exactStepper->RightHandSide(yInX, dydxRef);
            myStepper->RightHandSide(yIn, dydx);
            
        }
        
        cout.unsetf(ios_base::scientific);
        
        cout<<"\n# \t\t\t  -------| RESULTS |--------\n";
        cout<<"\n# \tStep_len = "<<step_len<<" mm "
        <<"\t theta_step = pi/"<<pi/theta;
        
        cout.setf(ios_base::scientific);
        cout.precision(2);
        
        cout<<"\n# The Max. error in each coordinate : \n#\t";
        for(int i=0;i <6; i++)
            cout<<maxErr[i]<<"  ";
        cout<<"\n The Average error in each coordinate : \n#\t";
        for(int i=0;i <6; i++)
            cout<<avgErr[i]/no_of_steps<<"  ";
        cout<<"\n\n#\t\t\t\t\t\t------------- End of output ----------------- \n\n";
        
    }
    
};




#endif /* defined(__Geant4__RKTest__) */


