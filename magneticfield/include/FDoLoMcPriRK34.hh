//  FDoLoMcPri4(3) header by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 7 July 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  FDoLoMcPriRK34.hh
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 7 July 2015


#ifndef FDoLo_McPri_34
#define FDoLo_McPri_34

#include "FSALMagIntegratorStepper.hh"

class FDoLoMcPriRK34 : public FSALMagIntegratorStepper
{
public:
    //Constructor using Equation
	FDoLoMcPriRK34(G4EquationOfMotion *EqRhs,
					 G4int numberOfVariables = 6,
					 G4bool primary =  true);
	~FDoLoMcPriRK34();

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
    G4int IntegratorOrder() const {return 3; }
    
private :
    
    FDoLoMcPriRK34(const FDoLoMcPriRK34&);
        //Copy constructor kept private
    FDoLoMcPriRK34& operator=(const FDoLoMcPriRK34&); 
        //assignment operator overloaded
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    G4double *pseudoDydx_for_DistChord ;
    
    FDoLoMcPriRK34* fAuxStepper;
};







#endif /* defined(__Geant4__FDoLoMcPriRK34__) */
