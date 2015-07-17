//  Bogacki-Shampine - 4 - 3(2) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 20 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  FBogackiShapine23.cc
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 20 May 2015
//
//
///////////////////////////////////////////////////////////////////////////////


#ifndef FBOGACKI_SHAMPINE23_H
#define FBOGACKI_SHAMPINE23_H

#include "FSALMagIntegratorStepper.hh"

class FBogackiShampine23 : public FSALMagIntegratorStepper{


 public:
 	//constructor
 	FBogackiShampine23( G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary= true ) ;

 	//destructor
 	~FBogackiShampine23() ;

 	//Stepper
 	 void Stepper(const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[],
                      	G4double nextDydx[]) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 2; }
     G4bool isFSAL() const{ return true; }
	G4double *getLastDydx();
    
	FBogackiShampine23(const FBogackiShampine23&);
   FBogackiShampine23& operator=(const FBogackiShampine23&);


 private:

	   G4double *ak2, *ak3, *ak4, *yTemp, *yIn;
      // for storing intermediate 'k' values in stepper
    G4double *pseudoDydx_for_DistChord;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations

    FBogackiShampine23* fAuxStepper;

//	G4int No_of_vars;
//	G4double hinit, tinit, tmax, *yinit;
//	double hmax, hmin, safe_const, err0, Step_factor;
//	void (*derivs)(double, double *, double *);


};

#endif /* FBogackiShampine23 */





























