#include "RKTest.hh"
//#include "RKTestDriver.hh"
//#include "src/RKTestDriver.cc"
#include "G4CashKarpRKF45.hh"
#include "VernerRK78.hh"

using namespace std;
int main(){

    RKTest myTest;
        int columns[] = {1,0,0,0,0,0};
//        myTest.testSteppersFixedUMF(columns);
    //    myTest.testPerformance("vr78");
    
    //    G4CashKarpRKF45 *theStpr(0);
    //    myTest.testAnyG4Stepper<G4CashKarpRKF45>("umf");
//        myTest.testSteppersFixedUMF<G4CashKarpRKF45>(columns);
        myTest.testStepperInterpolant<VernerRK78, G4CashKarpRKF45>(columns);
//    myTest.testG4Stepper<VernerRK78>("umf");
    //    cout<<"\n Hello";
    return 0;
}
