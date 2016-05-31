//
//  TsitourasRK45.hh
//  Geant4
//
//  Created by hackabot on 11/06/15.
//
//

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
    
    void SetupInterpolate( const G4double yInput[],
                              const G4double dydx[],
                              const G4double Step );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate( const G4double yInput[],
                         const G4double dydx[],
                         const G4double Step,
                         G4double yOut[],
                         G4double tau );
    
    
    void interpolate( const G4double yInput[],
                      const G4double dydx[],
                            G4double yOut[],
                            G4double Step,
                            G4double tau);

    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 4; }
    
//    G4double *getLastDydx();
    
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
