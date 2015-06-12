//  Verner - 9 - 6(5) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
//  Sponsored by Google in Google Summer of Code 2015.
// 
//  First version:  9 June 2015
//
//  This code is made available subject to the Geant4 license, a copy of
//  which is available at
//  http://geant4.org/license
//  
//  History
// -----------------------------
//  Created by Somnath on 9 June 2015.
//
//
/////////////////////////////////////////////////////////////////////////////


#ifndef VERNER_RK56_H
#define VERNER_RK56_H

#include "G4MagIntegratorStepper.hh"

class VernerRK56 : public G4MagIntegratorStepper{
    
    
public:
    //constructor
    VernerRK56( G4EquationOfMotion *EqRhs,
               G4int numberOfVariables = 6,
               G4bool primary= true ) ;
    
    //destructor
    ~VernerRK56() ;
    
    //Stepper
    void Stepper( const G4double y[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[] ) ;
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 5; }
     G4bool isFSAL() const{ return false; }
    
    VernerRK56(const VernerRK56&);
    VernerRK56& operator=(const VernerRK56&);
    
    
private:
    
	   G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *ak9,       
       // for storing intermediate 'k' values in stepper
    *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    VernerRK56* fAuxStepper;
    
    //	G4int No_of_vars;
    //	G4double hinit, tinit, tmax, *yinit;
    //	double hmax, hmin, safe_const, err0, Step_factor;
    //	void (*derivs)(double, double *, double *);
    
    
};

#endif /* VernerRK56 */

