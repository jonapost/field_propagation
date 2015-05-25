/*
 * BogackiShampine23.hh
 *
 *  Created on: 20-May-2015
 *      Author: hackabot
 */

#ifndef BOGACKI_SHAMPINE23_H
#define BOGACKI_SHAMPINE23_H

#include "G4MagIntegratorStepper.hh"

class BogackiShampine23 : public G4MagIntegratorStepper{


 public:
 	//constructor using equation
 	BogackiShampine23( G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary= true ) ;

 	//destructor
 	~BogackiShampine23() ;

 	//The Stepper
 	 void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 2; }

	BogackiShampine23(const BogackiShampine23&);
   BogackiShampine23& operator=(const BogackiShampine23&);


 private:

	   G4double *ak2, *ak3, *ak4, *yTemp, *yIn;
      // for storing intermediate 'k' values in stepper
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations

    BogackiShampine23* fAuxStepper;


};

#endif /* BogackiShampine23 */





























