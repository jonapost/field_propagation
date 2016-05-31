//
//  FTsitourasRK45.hh
//  Geant4
//
//  Created by hackabot on 1 August 2015.
//
//

#ifndef FTsitouras_RK45
#define FTsitouras_RK45

#include "FSALMagIntegratorStepper.hh"

class FTsitourasRK45 : public FSALMagIntegratorStepper
{
public:
    FTsitourasRK45(G4EquationOfMotion *EqRhs,
                  G4int numberOfVariables = 6,
                  G4bool primary =  true);
    ~FTsitourasRK45();
    
    void Stepper( const G4double y[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[],
                 G4double nextDydx[] ) ;
    
    void interpolate( const G4double yInput[],
                     const G4double dydx[],
                     G4double yOut[],
                     G4double Step,
                     G4double tau);
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 4; }
    G4bool isFSAL() const{ return false; }
    
    //    G4double *getLastDydx();
    
    private :
    
    FTsitourasRK45(const FTsitourasRK45&);
    FTsitourasRK45& operator=(const FTsitourasRK45&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *yTemp, *yIn;
    
    G4double *fake_Next_Dydx;
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    FTsitourasRK45* fAuxStepper;
};



#endif /* defined(__Geant4__FTsitourasRK45__) */
