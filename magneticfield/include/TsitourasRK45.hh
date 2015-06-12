//  Tsitouras- 7 - 5(4) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
//  Sponsored by Google in Google Summer of Code 2015.
// 
//  First version:  11 June 2015
//
//  This code is made available subject to the Geant4 license, a copy of
//  which is available at
//  http://geant4.org/license
//  
//  History
// -----------------------------
//  Created by Somnath on 11 June 2015.
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef Tsitouras_RK45
#define Tsitouras_RK45

#include "G4MagIntegratorStepper.hh"

class TsitourasRK45 : public G4MagIntegratorStepper
{
public:
	TsitourasRK45(G4EquationOfMotion *EqRhs,
					 G4int numberOfVariables = 6,
					 G4bool primary =  true);
	~TsitourasRK45();

	void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 4; }
    
private :
    
    TsitourasRK45(const TsitourasRK45&);
    TsitourasRK45& operator=(const TsitourasRK45&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    TsitourasRK45* fAuxStepper;
};



#endif /* defined(__Geant4__TsitourasRK45__) */
