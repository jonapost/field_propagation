//  Bogacki-Shampine - 8 - 5(4) implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 26 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  DormandPrince745.cc
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 26 May 2015
//

#ifndef Bogacki_Shampine_45
#define Bogacki_Shampine_45

#include "G4MagIntegratorStepper.hh"

class BogackiShampine45 : public G4MagIntegratorStepper
{
public:
    BogackiShampine45(G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary =  true);
    ~BogackiShampine45();
    
    void Stepper( const G4double y[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[] ) ;
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 4; }
    
    private :
    
    BogackiShampine45(const BogackiShampine45&);
    BogackiShampine45& operator=(const BogackiShampine45&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    BogackiShampine45* fAuxStepper;
};







#endif /* defined(__Geant4__BogackiShampine45__) */
