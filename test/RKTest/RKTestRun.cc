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
//Old G4 Steppers
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"


using namespace std;
int main(int argc, char* argv[]){
    
    RKTest myTest;
        int columns[] = {0,0,0,0,0,0};
    
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

    
    cout<<" ::::: INTERPOLATION ACCURACY RECORD ::: ";
    cout<<"\n\n TESTING PERFORMANCE IN UMF \n";
    cout<< "\n    :::> Field_factor = "<<Field_factor<<endl;

    G4double div_list[] = {0.5,2.0/3.0,3.0/4.0,1.0,2.0,3.0,4.0,6.0,8.0,12.0 };

    int i=0;
    while(div_list[i] && i<50){
        cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
        <<"Testing VernerRK78 ...\n" ;
        myTest.testStepperInterpolant<VernerRK78, VernerRK78>(columns,div_list[i], 4.0*pi) ;
        i++;
    }

//    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
//    <<"Testing VernerRK67  ...\n" ;
//    myTest.testStepperInterpolant<VernerRK67, VernerRK78>(columns); ;
//    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
//    <<"Testing VernerRK78  ...\n" ;
//    myTest.testStepperInterpolant<VernerRK78, VernerRK78>(columns); ;
//    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
//    <<"Testing BogackiShampine45  ...\n" ;
//    myTest.testStepperInterpolant<BogackiShampine45, VernerRK78>(columns); ;
//    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
//    <<"Testing DormandPrince745  ...\n" ;
//    myTest.testStepperInterpolant<DormandPrince745, VernerRK78>(columns); ;
////    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
////    <<"Testing BogackiShampine23  ...\n" ;
////    myTest.testStepperInterpolant<BogackiShampine23, VernerRK78>(columns); ;
//    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
//    <<"Testing DormandPrinceRK56  ...\n" ;
//    myTest.testStepperInterpolant<DormandPrinceRK56, VernerRK78>(columns); ;
////    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
////    <<"Testing DormandPrinceRK78  ...\n" ;
////    myTest.testStepperInterpolant<DormandPrinceRK78, VernerRK78>(columns); ;
//    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
//    <<"Testing DoLoMcPriRK34  ...\n" ;
//    myTest.testStepperInterpolant<DoLoMcPriRK34, VernerRK78>(columns); ;
//    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
//    <<"Testing TsitourasRK45  ...\n" ;
//    myTest.testStepperInterpolant<TsitourasRK45, VernerRK78>(columns);
////    cout<<"\n\n :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: \n"
////    <<"Testing G4CashKarpRKF45  ...\n" ;
////    myTest.testStepperInterpolant<G4CashKarpRKF45, VernerRK78>(columns);

    return 0;

}
