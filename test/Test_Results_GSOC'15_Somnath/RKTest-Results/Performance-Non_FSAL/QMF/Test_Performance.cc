#include "RKTest.hh"
//#include "RKTestDriver.hh"
//#include "src/RKTestDriver.cc"
#include "VernerRK56.hh"
#include "VernerRK67.hh"
#include "VernerRK78.hh"
#include "BogackiShampine45.hh"
#include "DormandPrince745.hh"
#include "BogackiShampine23.hh"
#include "DormandPrinceRK56.hh"
#include "DormandPrinceRK78.hh"
#include "DoLoMcPriRK34.hh"
#include "TsitourasRK45.hh"

using namespace std;
int main(int argc, char* argv[]){

    RKTest myTest;
//    int columns[] = {1,1,0,0,0,0};

    G4double Field_factor = 1.0;
    
    if(argc > 1)
        Field_factor = atof(argv[1]);
    
    
    //    G4CashKarpRKF45 *theStpr(0);
//        myTest.testG4Stepper<G4CashKarpRKF45>("umf",10);
//    myTest.testFSALStepper<FBogackiShampine45>("umf", 10);
//        myTest.testSteppersFixedUMF<G4CashKarpRKF45>(columns);
//        myTest.testStepperInterpolant<TsitourasRK45,
//        					BogackiShampine45>(columns);
//    myTest.testG4Stepper<VernerRK78>("umf");
    
    

//    myTest.testSteppersFixedUMF<<#class STEPPER#>>(<#int *columns#>)

    cout<<" ::::: NON-FSAL STEPPERS PERFORMANCE RECORD ::: ";
    cout<<"\n\n TESTING PERFORMANCE IN QMF \n";
    cout<< "\n    :::> Field_factor = "<<Field_factor<<endl;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing VernerRK56 ...\n" ;
    myTest.testG4Stepper<VernerRK56>("qmf", Field_factor) ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing VernerRK67  ...\n" ;
    myTest.testG4Stepper<VernerRK67>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing VernerRK78  ...\n" ;
    myTest.testG4Stepper<VernerRK78>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing BogackiShampine45  ...\n" ;
    myTest.testG4Stepper<BogackiShampine45>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing DormandPrince745  ...\n" ;
    myTest.testG4Stepper<DormandPrince745>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing BogackiShampine23  ...\n" ;
    myTest.testG4Stepper<BogackiShampine23>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing DormandPrinceRK56  ...\n" ;
    myTest.testG4Stepper<DormandPrinceRK56>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing DormandPrinceRK78  ...\n" ;
    myTest.testG4Stepper<DormandPrinceRK78>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing DoLoMcPriRK34  ...\n" ;
    myTest.testG4Stepper<DoLoMcPriRK34>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing TsitourasRK45  ...\n" ;
    myTest.testG4Stepper<TsitourasRK45>("qmf", Field_factor); ;
    return 0;
}
