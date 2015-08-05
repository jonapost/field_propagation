


#include "G4ExactHelixStepper.hh"	//For exact UMF stepper
#include "RKTestDriver.hh"
#include "FSALMagIntegratorDriver.hh"

template<class STEPPER>
G4MagInt_Driver *
SetDriver(/*std::string my_stepper_code,*/ G4Mag_UsualEqRhs *fEquation /*G4FieldTrack myTrack*/ ){
    
    G4double hminimum = 1.e-12;
    STEPPER *theStepper = new STEPPER(fEquation);
    G4MagInt_Driver *theDriver = new G4MagInt_Driver(hminimum, theStepper);
    return theDriver;
    
}

template <class STEPPER>
FSALMagInt_Driver *SetFSALDriver(G4Mag_UsualEqRhs *fEquation){
    G4double hminimum = 1.e-12;
    STEPPER *theStepper = new STEPPER(fEquation);
    
    
    FSALMagInt_Driver *theDriver = new FSALMagInt_Driver(hminimum, theStepper);
    
//    FSALMagInt_Driver *theDriver = new FSALMagInt_Driver(hminimum,theStepper);//, 0, 0);
    return theDriver;
}