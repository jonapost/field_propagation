//  DoLoMcPri4(3) header by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 7 July 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  DoLoMcPriRK34.hh
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 7 July 2015


#ifndef DoLo_McPri_34
#define DoLo_McPri_34

#include "G4MagIntegratorStepper.hh"

class DoLoMcPriRK34 : public G4MagIntegratorStepper
{
public:
    //Constructor using Equation
	DoLoMcPriRK34(G4EquationOfMotion *EqRhs,
					 G4int numberOfVariables = 6,
					 G4bool primary =  true);
	~DoLoMcPriRK34();

    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;
    
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
    
    void interpolate(const G4double yInput[],
                     const G4double dydx[],
                           G4double yOut[],
                           G4double Step,
                           G4double tau ) ;


    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 3; }
    
private :
    
    DoLoMcPriRK34(const DoLoMcPriRK34&);
        //Copy constructor kept private
    DoLoMcPriRK34& operator=(const DoLoMcPriRK34&); 
        //assignment operator overloaded
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    DoLoMcPriRK34* fAuxStepper;
};







#endif /* defined(__Geant4__DoLoMcPriRK34__) */
