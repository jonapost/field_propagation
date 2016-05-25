//  Bogacki-Shampine - 8 - 5(4) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 26 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  FBogackiShampine45.hh
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 26 May 2015
//
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FBogacki_Shampine_45
#define FBogacki_Shampine_45

#include "FSALMagIntegratorStepper.hh"

class FBogackiShampine45 : public FSALMagIntegratorStepper
{
public:
    FBogackiShampine45(G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary =  true);
    ~FBogackiShampine45();
    
    void Stepper( const G4double y[],
                  const G4double dydx[],
              	   		G4double h,
                	 	G4double yout[],
                 		G4double yerr[],
                 		G4double nextDydx[]) ;
    
    void interpolate( const G4double yInput[],
                      const G4double dydx[],
                            G4double yOut[],
                            G4double Step,
                            G4double tau ) ;
    
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 4; }
    
    private :
    
    FBogackiShampine45(const FBogackiShampine45&);
    FBogackiShampine45& operator=(const FBogackiShampine45&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8,
    *ak9, *ak10, *ak11,
    *yTemp, *yIn;
    
    G4double *pseudoDydx_for_DistChord;
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    FBogackiShampine45* fAuxStepper;
};

#endif /* defined(__Geant4__FBogackiShampine45__) */
