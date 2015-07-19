//
//  RKTest.h
//  Geant4
//
//  Created by Somnath on 14/07/15.
//
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
#include <stdio.h>
using namespace std;

class RKTest{
private:
    //Members - all needed for the program Run
    G4UniformMagField uField ;
    G4QuadrupoleMagField QField ;
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
    
//    void testSteppersFixed(string field_code, string stepper_code, G4double step_len);
    
    template < class STEPPER >
    void testSteppersFixedUMF(int columns[6],
                              /*string stepper_code = "ck45",*/
                              G4double step_len = 25.0*CLHEP::mm,
                              int no_of_steps = 100);
    
//    void testSteppersFixedQMF(string stepper_code = "ck45", G4double step_len = 25.0*CLHEP::mm);

    template<class STEPPER>
    void testG4Stepper(string field_code = "umh");
    
    
    template<class STEPPER>
    void testFSALStepper(string field_code = "umf");
    
    template<class STEPPER, class REF_STEPPER>
    void testStepperInterpolant(int columns[6], string field_code = "umf", G4double step_len_pi_divisor = 6.0, G4double maxAngle = 2.0*CLHEP::pi);
    
    static void print3(int columns[], int ifHeader , G4double yOut[], G4double yErr[], G4double yOutX[]);
//    void print3Data(G4double yOut[], yErr[], yOutX);
    
    
};



#endif /* defined(__Geant4__RKTest__) */
