//
//  RKTest.cpp
//  Geant4
//
//  Created by Somnath on 14/07/15.
//
//

#include "RKTest.hh"
//#include "RKTestDriver.cc"
//#include <iostream>
#include <iomanip>
using namespace CLHEP;
using namespace std;

RKTest::RKTest() :

tTrack(G4FieldTrack(CreateTrack())){

}

RKTest::~RKTest(){
    //Do nothing
}



////    tDriver = SetDriver(stepper_code, QField, tTrack);
//    //What Exact/Solution stepper are we gonna use ??
//}


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



void RKTest::Reset(){
//    uField = SetupUMF();
//    QField = SetupQMF();

    if(MField)
        delete MField;
    tTrack = CreateTrack();
    if(tDriver)
        delete tDriver;
    if(fEqRhs)
        delete fEqRhs;
}

