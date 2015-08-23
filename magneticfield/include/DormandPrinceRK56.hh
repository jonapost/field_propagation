//  Dormand-Prince RK 6(5) non-FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 26 June 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  DormandPrince745.cc
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 26 June 2015
//
//
///////////////////////////////////////////////////////////////////////////////

#ifndef DORMAND_PRINCE_RK56_H
#define DORMAND_PRINCE_RK56_H

#include "G4MagIntegratorStepper.hh"

class DormandPrinceRK56 : public G4MagIntegratorStepper{
    
    
public:
    //constructor
    DormandPrinceRK56( G4EquationOfMotion *EqRhs,
               G4int numberOfVariables = 6,
               G4bool primary= true ) ;
    
    //destructor
    ~DormandPrinceRK56() ;
    
    //Stepper
    void Stepper( const G4double y[],
                  const G4double dydx[],
                 		G4double h,
                 		G4double yout[],
                 		G4double yerr[] ) ;
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 5; }
    


    

    
    DormandPrinceRK56(const DormandPrinceRK56&);
    DormandPrinceRK56& operator=(const DormandPrinceRK56&);
    
    
    //For Preparing the Interpolant and calculating the extra stages
    void SetupInterpolate_low( const G4double yInput[],
                              const G4double dydx[],
                              const G4double Step );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate_low( const G4double yInput[],
                         const G4double dydx[],
                         const G4double Step,
                         G4double yOut[],
                         G4double tau );
    
    inline void SetupInterpolate( const G4double yInput[],
                                 const G4double dydx[],
                                 const G4double Step ){
        SetupInterpolate_low( yInput, dydx, Step);
    }
    
    //For calculating the output at the tau fraction of Step
    inline void Interpolate( const G4double yInput[],
                            const G4double dydx[],
                            const G4double Step,
                            G4double yOut[],
                            G4double tau ){
        Interpolate_low( yInput, dydx, Step, yOut, tau);
    }
    
    void SetupInterpolate_high( const G4double yInput[],
                               const G4double dydx[],
                               const G4double Step );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate_high( const G4double yInput[],
                          const G4double dydx[],
                          const G4double Step,
                          G4double yOut[],
                          G4double tau );
    
    
    
private:
    
	   G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *ak9,       // for storing intermediate 'k' values in stepper
    *ak10_low, *ak10, *ak11, * ak12, //For the additional stages of Interpolant
    *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    DormandPrinceRK56* fAuxStepper;
    
    
    
};

#endif /* DormandPrinceRK56 */

