//RKTestDriver.hh
//Created by Somnath on 14 July 2015
//Sponsored by Google Summer of Code

#include "FSALMagIntegratorDriver.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "FSALMagIntegratorStepper.hh"
#include <string>

//Include all the steppers
#include "G4ExplicitEuler.hh"
#include "VernerRK67.hh"
#include "VernerRK56.hh"
#include "VernerRK78.hh"
#include "G4CashKarpRKF45.hh"
#include "G4ClassicalRK4.hh"
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
//template<typename MAG_FIELD>
G4MagInt_Driver *
SetDriver(std::string my_stepper_code, G4Mag_UsualEqRhs *fEquation, G4FieldTrack myTrack );
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
