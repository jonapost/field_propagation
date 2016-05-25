/*
 * VernerRK67.hh
 *
 *  Created on: 8-Jun-2015
 *      Author: hackabot
 */

#ifndef VERNER_RK67_H
#define VERNER_RK67_H

#include "G4MagIntegratorStepper.hh"

class VernerRK67 : public G4MagIntegratorStepper{


 public:
 	//constructor
 	VernerRK67( G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary= true ) ;

 	//destructor
 	~VernerRK67() ;

 	//Stepper
 	 void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    
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
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 6; }

	VernerRK67(const VernerRK67&);
   VernerRK67& operator=(const VernerRK67&);


 private:

	   G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *ak9, *ak10,       // for storing intermediate 'k' values in stepper
    
    *ak11, *ak12, *ak13, *ak14, *ak15, *ak16, // for use in the interpolants
     *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations

    VernerRK67* fAuxStepper;

//	G4int No_of_vars;
//	G4double hinit, tinit, tmax, *yinit;
//	double hmax, hmin, safe_const, err0, Step_factor;
//	void (*derivs)(double, double *, double *);


};

#endif /* VernerRK67 */

