
//
//  DormandPrinceRK78.hh
//  Geant4
//
//  Created by Somnath on 30/06/15.
//
//

#ifndef Dormand_Prince_RK78
#define Dormand_Prince_RK78

#include "G4MagIntegratorStepper.hh"


class DormandPrinceRK78 : public G4MagIntegratorStepper
{
public:
    DormandPrinceRK78(G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary =  true);
    ~DormandPrinceRK78();
    
    void Stepper( const G4double y[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[]) ;
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 7; }
    
    private :
    
    DormandPrinceRK78(const DormandPrinceRK78&);
    DormandPrinceRK78& operator=(const DormandPrinceRK78&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8,
    *ak9, *ak10, *ak11, *ak12, *ak13,
    *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    DormandPrinceRK78* fAuxStepper;
};







#endif /* defined(__Geant4__DormandPrinceRK78__) */
