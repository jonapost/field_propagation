////Include all the steppers
//#include "G4ExplicitEuler.hh"
//#include "VernerRK67.hh"
//#include "VernerRK56.hh"
//#include "VernerRK78.hh"
//#include "G4CashKarpRKF45.hh"
//#include "G4ClassicalRK4.hh"
//#include "G4SimpleHeum.hh"
//#include "DormandPrince745.hh"
//#include "BogackiShampine45.hh"
//#include "G4ExactHelixStepper.hh"
//#include "BogackiShampine23.hh"
////Include class header
//#include "RKTestDriver.hh"
////template <class INT_DRIVER, class MAG_FIELD, class STEPPER >
//////class RKTestDriver
//////
//////void RKTestDriver<class INT_DRIVER,  INT_DRIVER, MAG_FIELD, INT_STEPPER > ::
////template<typename MAG_FIELD>
////G4MagInt_Driver *
////SetDriver(std::string my_stepper_code, MAG_FIELD myField, G4FieldTrack myTrack ){
////    
////    //    strcpy(my_stepper_code, my_stepper_code);
////    
////    G4double hminimum = 1.e-12;
////    G4MagIntegratorStepper *theStepper;
////    std::string stepper_name;
////    
////    G4double momentum = myTrack.GetMomentum().mag();
////    G4ChargeState chargestate = *myTrack.GetChargeState();
////    G4double mass = myTrack.GetRestMass();
////    
////    
////    G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);
////    fEquation->SetChargeMomentumMass(chargestate, momentum, mass);
////    
////    //Matching Stepper code to assign the stepper
////    if( my_stepper_code == "ck45") {
////        theStepper = new G4CashKarpRKF45(fEquation);
////        stepper_name = "G4CashKarpRKF45";
////    }
////    //        if(strcmp( my_stepper_code,  0)) {
////    //            stages = G4ExactHelix_stages - 1;
////    //            stepper_name = "Exact";
////    //            }
////    if( my_stepper_code ==  "bs23" ) { theStepper = new BogackiShampine23(fEquation);
////        stepper_name = "BogackiShampine23";
////    }
////    if( my_stepper_code ==  "dp45" ) { theStepper = new DormandPrince745(fEquation);
////        stepper_name = "DormandPrince745";
////    }
////    if( my_stepper_code ==  "bs45" ) { theStepper = new BogackiShampine45(fEquation);
////        stepper_name = "BogackiShampine45";
////    }
////    if( my_stepper_code ==  "clc4" ) { theStepper = new G4ClassicalRK4(fEquation);
////        stepper_name = "G4ClassicalRK4";
////    }
////    //        if( my_stepper_code ==  6) { theStepper = new G4SimpleHeum(fEquation);
////    //            stages = G4SimpleHeum_stages - 1;
////    //            stepper_name = "G4SimpleHeum";
////    //            }
////    //        if( my_stepper_code ==  7) { theStepper = new G4ExplicitEuler(fEquation);
////    //            stages = G4ExplicitEuler_stages - 1;
////    //            stepper_name = "G4ExplicitEuler";
////    //            }
////    if( my_stepper_code ==  "vr67" ) { theStepper = new VernerRK67(fEquation);
////        stepper_name = "VernerRK67";
////    }
////    if( my_stepper_code ==  "vr56") { theStepper = new VernerRK56(fEquation);
////        stepper_name = "VernerRK56";
////    }
////    if( my_stepper_code ==  "vr78") { theStepper = new VernerRK78(fEquation);
////        stepper_name = "VernerRK78";
////    }
////    
////    G4MagInt_Driver *theDriver = new G4MagInt_Driver(hminimum, theStepper);
////    return theDriver;
////    
////}
////
////
//////template<class INT_DRIVER, class MAG_FIELD, class INT_STEPPER >
//////void RKTestDriver<INT_DRIVER, MAG_FIELD, INT_STEPPER >
//////::ResetDriver(){
//////    SetDriver(my_stepper_code, fEquation->GetFieldObj());
//////}
////
////
//////template<class INT_DRIVER, class MAG_FIELD, class INT_STEPPER >
//////RKTestDriver<INT_DRIVER, MAG_FIELD, INT_STEPPER > ::RKTestDriver(){
//////    //Do nothing
//////}
////
//////template<class INT_DRIVER, class MAG_FIELD, class INT_STEPPER >
//////void RKTestDriver<INT_DRIVER, MAG_FIELD, INT_STEPPER >
//////::
//////template<class INT_DRIVER, class MAG_FIELD, class INT_STEPPER >
//////void RKTestDriver<INT_DRIVER, MAG_FIELD, INT_STEPPER >
//////
//////template<class INT_DRIVER, class MAG_FIELD, class INT_STEPPER >
//////void RKTestDriver<INT_DRIVER, MAG_FIELD, INT_STEPPER >
//////
//////
//////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
//
//#include "RKTestDriver.hh"
//
//template<class STEPPER>
//G4MagInt_Driver *
//SetDriver(/*std::string my_stepper_code,*/ G4Mag_UsualEqRhs *fEquation, G4FieldTrack myTrack ){
//    
//    //    strcpy(my_stepper_code, my_stepper_code);
//    
//    G4double hminimum = 1.e-12;
//    STEPPER *theStepper = new STEPPER(fEquation);
////    std::string stepper_name;
//    
//    //    G4double momentum = myTrack.GetMomentum().mag();
//    //    G4ChargeState chargestate = *myTrack.GetChargeState();
//    //    G4double mass = myTrack.GetRestMass();
//    //
//    //
//    //    G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);
//    //    fEquation->SetChargeMomentumMass(chargestate, momentum, mass);
//    //    G4Mag_UsualEqRhs *fEquation = SetMagEqRhs(myField, myTrack);
//    //Matching Stepper code to assign the stepper
////    if( my_stepper_code == "ck45") {
////        theStepper = new G4CashKarpRKF45(fEquation);
////        stepper_name = "G4CashKarpRKF45";
////    }
////    //        if(strcmp( my_stepper_code,  0)) {
////    //            stages = G4ExactHelix_stages - 1;
////    //            stepper_name = "Exact";
////    //            }
////    if( my_stepper_code ==  "bs23" ) { theStepper = new BogackiShampine23(fEquation);
////        stepper_name = "BogackiShampine23";
////    }
////    if( my_stepper_code ==  "dp45" ) { theStepper = new DormandPrince745(fEquation);
////        stepper_name = "DormandPrince745";
////    }
////    if( my_stepper_code ==  "bs45" ) { theStepper = new BogackiShampine45(fEquation);
////        stepper_name = "BogackiShampine45";
////    }
////    if( my_stepper_code ==  "clc4" ) { theStepper = new G4ClassicalRK4(fEquation);
////        stepper_name = "G4ClassicalRK4";
////    }
////    //        if( my_stepper_code ==  6) { theStepper = new G4SimpleHeum(fEquation);
////    //            stages = G4SimpleHeum_stages - 1;
////    //            stepper_name = "G4SimpleHeum";
////    //            }
////    //        if( my_stepper_code ==  7) { theStepper = new G4ExplicitEuler(fEquation);
////    //            stages = G4ExplicitEuler_stages - 1;
////    //            stepper_name = "G4ExplicitEuler";
////    //            }
////    if( my_stepper_code ==  "vr67" ) { theStepper = new VernerRK67(fEquation);
////        stepper_name = "VernerRK67";
////    }
////    if( my_stepper_code ==  "vr56") { theStepper = new VernerRK56(fEquation);
////        stepper_name = "VernerRK56";
////    }
////    if( my_stepper_code ==  "vr78") { theStepper = new VernerRK78(fEquation);
////        stepper_name = "VernerRK78";
////    }
//    
//    G4MagInt_Driver *theDriver = new G4MagInt_Driver(hminimum, theStepper);
//    return theDriver;
//    
//}





//Include all the steppers
#include "G4ExplicitEuler.hh"
#include "VernerRK67.hh"
#include "VernerRK56.hh"
#include "VernerRK78.hh"
#include "G4CashKarpRKF45.hh"
//#include "G4ClassicalRK4.hh"
#include "G4SimpleHeum.hh"
#include "DormandPrince745.hh"
#include "BogackiShampine45.hh"
#include "G4ExactHelixStepper.hh"
#include "BogackiShampine23.hh"
//Include class header
//#include "RKTestDriver.hh"
//template <class INT_DRIVER, class MAG_FIELD, class STEPPER >
//class RKTestDriver
//
//void RKTestDriver<class INT_DRIVER,  INT_DRIVER, MAG_FIELD, INT_STEPPER > ::
//template<class STEPPER>
//G4MagInt_Driver *
//SetDriver(/*std::string my_stepper_code ,*/ G4Mag_UsualEqRhs *fEquation, G4FieldTrack myTrack );
//template<typename MAG_FIELD>
//G4Mag_UsualEqRhs *SetMagEqRhs(MAG_FIELD myField, G4FieldTrack myTrack ){
//    G4double momentum = myTrack.GetMomentum().mag();
//    G4ChargeState chargestate = *myTrack.GetChargeState();
//    G4double mass = myTrack.GetRestMass();
//
//
//    G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&myField);
//    fEquation->SetChargeMomentumMass(chargestate, momentum, mass);
//    return fEquation;
//}


#include "RKTestDriver.hh"

template<class STEPPER>
G4MagInt_Driver *
SetDriver(/*std::string my_stepper_code,*/ G4Mag_UsualEqRhs *fEquation, G4FieldTrack myTrack ){
    
    G4double hminimum = 1.e-12;
    STEPPER *theStepper = new STEPPER(fEquation);
    G4MagInt_Driver *theDriver = new G4MagInt_Driver(hminimum, theStepper);
    return theDriver;
    
}
