/*
 * BogackiShampine23.hh
 *
 *  Created on: 20-May-2015
 *      Author: hackabot
 */

#ifndef MODIFIED_BOGACKI_SHAMPINE23_H
#define MODIFIED_BOGACKI_SHAMPINE23_H

#include "G4MagIntegratorStepper.hh"

class BogackiShampine23modified : public G4MagIntegratorStepper{


 public:
 	//constructor using equation
 	BogackiShampine23modified( G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary= true ) ;

 	//destructor
 	~BogackiShampine23modified() ;

 	//The Stepper
 	 void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 2; }

	BogackiShampine23modified(const BogackiShampine23modified&);
   BogackiShampine23modified& operator=(const BogackiShampine23modified&);
   
   void modifiedRightHandSide(const G4double yIn[], G4double dydx[]);


 private:

	   G4double *ak2, *ak3, *ak4, *yTemp, *yIn;
      // for storing intermediate 'k' values in stepper
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations

    BogackiShampine23modified* fAuxStepper;


};

#endif /* BogackiShampine23modified */





























