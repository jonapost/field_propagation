//  Verner's RK 6(5) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 9 June 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  FVernerRK56.hh
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 9 June 2015
//
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FVERNER_RK56_H
#define FVERNER_RK56_H

#include "G4VFSALIntegrationStepper.hh"

class FVernerRK56 : public G4VFSALIntegrationStepper{
    
    
public:
    //constructor
    FVernerRK56( G4EquationOfMotion *EqRhs,
               G4int numberOfVariables = 6,
               G4bool primary= true ) ;
    
    //destructor
    ~FVernerRK56() ;
    
    //Stepper
    void Stepper( const G4double yIn[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[],
                 G4double nextDydx[] ) ;
    
    void interpolate( const G4double yInput[],
                     const G4double dydx[],
                     G4double yOut[],
                     G4double Step,
                     G4double tau
                     ) ;
    
    void interpolate6( const G4double yInput[],
                      const G4double dydx[],
                      G4double yOut[],
                      G4double Step,
                      G4double tau
                      ) ;
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 5; }

    
    FVernerRK56(const FVernerRK56&);
    FVernerRK56& operator=(const FVernerRK56&);
    
    
private:
    
	   G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *ak9,       // for storing intermediate 'k' values in stepper
    *ak10, *ak11, *ak12, //for storing extra stages in interpolate 
    *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    FVernerRK56* fAuxStepper;
    G4double *pseudoDydx_for_DistChord ;

    
    //	G4int No_of_vars;
    //	G4double hinit, tinit, tmax, *yinit;
    //	double hmax, hmin, safe_const, err0, Step_factor;
    //	void (*derivs)(double, double *, double *);
    
    
};

#endif /* FVernerRK56 */

