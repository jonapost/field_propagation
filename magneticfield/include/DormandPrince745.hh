//  DormandPrince7 - 5(4) implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 25 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  DormandPrince745.cc
//  Geant4
//
//  History
// ------------------------------------------
//  Created : 25 May 2015.             - Somnath
//
//  Added interpolate() :               - Somnath
//      29 June 2015


#ifndef Dormand_Prince_745
#define Dormand_Prince_745

#include "G4MagIntegratorStepper.hh"

class DormandPrince745 : public G4MagIntegratorStepper
{
public:
	DormandPrince745(G4EquationOfMotion *EqRhs,
					 G4int numberOfVariables = 6,
					 G4bool primary =  true);
	~DormandPrince745();

	void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;
    
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

    G4double  DistChord() const;
    G4double DistChord2() const;
    G4int IntegratorOrder() const {return 4; }
    
    //Newly constructed copy constructor
//    DormandPrince745(const DormandPrince745&);
    
private :
    
    DormandPrince745& operator=(const DormandPrince745&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7,
    *ak8, *ak9, 	//For additional stages in the interpolant
    *yTemp, *yIn;
    
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    DormandPrince745* fAuxStepper;
};







#endif /* defined(__Geant4__DormandPrince745__) */
