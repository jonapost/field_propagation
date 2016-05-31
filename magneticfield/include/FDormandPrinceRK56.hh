//  FDormand-Prince RK 6(5) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 26 June 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  FDormandPrince745.hh
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 26 June 2015
//
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FDORMAND_PRINCE_RK56_H
#define FDORMAND_PRINCE_RK56_H

#include "FSALMagIntegratorStepper.hh"

class FDormandPrinceRK56 : public FSALMagIntegratorStepper{
    
    
public:
    //constructor
    FDormandPrinceRK56( G4EquationOfMotion *EqRhs,
               G4int numberOfVariables = 6,
               G4bool primary= true ) ;
    
    //destructor
    ~FDormandPrinceRK56() ;
    
    //Stepper
    void Stepper( const G4double y[],
                  const G4double dydx[],
                 		G4double h,
                 		G4double yout[],
                 		G4double yerr[],
                        G4double nextDydx[] ) ;
    
    void interpolate( const G4double yInput[],
                      const G4double dydx[],
                      const G4double Step,
                     		G4double yOut[],
                     		G4double tau );
    
    //For Preparing the Interpolant and calculating the extra stages
    void SetupInterpolate( const G4double yInput[],
                           const G4double dydx[],
                           const G4double Step );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate( const G4double yInput[],
                      const G4double dydx[],
                      const G4double Step,
                     		G4double yOut[],
                     		G4double tau );
    

    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 5; }
//    G4bool isFSAL() const{ return false; }
    //    G4double *getLastDydx() {return 0;};
    
    FDormandPrinceRK56(const FDormandPrinceRK56&);
    FDormandPrinceRK56& operator=(const FDormandPrinceRK56&);
    
    
private:
    
	   G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *ak9,       // for storing intermediate 'k' values in stepper
    *ak10, *ak11, * ak12, //For the additional stages of Interpolant
    *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    
    G4double *pseudoDydx_for_DistChord;
    // for DistChord calculations
    
    FDormandPrinceRK56* fAuxStepper;
    
    //	G4int No_of_vars;
    //	G4double hinit, tinit, tmax, *yinit;
    //	double hmax, hmin, safe_const, err0, Step_factor;
    //	void (*derivs)(double, double *, double *);
    
    
};

#endif /* FDormandPrinceRK56 */

