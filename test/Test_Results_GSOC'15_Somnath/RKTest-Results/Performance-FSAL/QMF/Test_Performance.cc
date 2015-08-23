#include "RKTest.hh"
//#include "RKTestDriver.hh"
//#include "src/RKTestDriver.cc"
#include "FVernerRK56.hh"
#include "FBogackiShampine45.hh"
#include "FDormandPrince745.hh"
#include "FBogackiShampine23.hh"
#include "FDormandPrinceRK56.hh"
#include "FDoLoMcPriRK34.hh"
#include "FTsitourasRK45.hh"

using namespace std;
int main(int argc, char* argv[]){

    RKTest myTest;
//    int columns[] = {1,1,0,0,0,0};

    G4double Field_factor = 1.0;
    
    if(argc > 1)
        Field_factor = atof(argv[1]);
    
    
    //    G4CashKarpRKF45 *theStpr(0);
//        myTest.testFSALStepper<G4CashKarpRKF45>("umf",10);
//    myTest.testFSALStepper<FBogackiShampine45>("umf", 10);
//        myTest.testSteppersFixedUMF<G4CashKarpRKF45>(columns);
//        myTest.testStepperInterpolant<TsitourasRK45,
//                          BogackiShampine45>(columns);
//    myTest.testFSALStepper<VernerRK78>("umf");
    
    

//    myTest.testSteppersFixedUMF<<#class STEPPER#>>(<#int *columns#>)

    cout<<" ::::: NON-FSAL STEPPERS PERFORMANCE RECORD ::: ";
    cout<<"\n\n TESTING PERFORMANCE IN QMF \n";
    cout<< "\n    :::> Field_factor = "<<Field_factor<<endl;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing FVernerRK56 ...\n" ;
    myTest.testFSALStepper<FVernerRK56>("qmf", Field_factor) ;

    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing FBogackiShampine45  ...\n" ;
    myTest.testFSALStepper<FBogackiShampine45>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing FDormandPrince745  ...\n" ;
    myTest.testFSALStepper<FDormandPrince745>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing FBogackiShampine23  ...\n" ;
    myTest.testFSALStepper<FBogackiShampine23>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing FDormandPrinceRK56  ...\n" ;
    myTest.testFSALStepper<FDormandPrinceRK56>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing FDoLoMcPriRK34  ...\n" ;
    myTest.testFSALStepper<FDoLoMcPriRK34>("qmf", Field_factor); ;
    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
    <<"Testing FTsitourasRK45  ...\n" ;
    myTest.testFSALStepper<FTsitourasRK45>("qmf", Field_factor); ;
    return 0;
}
