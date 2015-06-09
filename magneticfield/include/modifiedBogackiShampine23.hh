/*
 * BogackiShampine23.hh
 *
 *  Created on: 20-May-2015
 *      Author: hackabot
 */

#ifndef MODIFIED_BOGACKI_SHAMPINE23_H
#define MODIFIED_BOGACKI_SHAMPINE23_H

#include "G4MagIntegratorStepper.hh"

class modifiedBogackiShampine23 : public G4MagIntegratorStepper{


 public:
 	//constructor using equation
 	modifiedBogackiShampine23( G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary= true ) ;

 	//destructor
 	~modifiedBogackiShampine23() ;

 	//The Stepper
 	 void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 2; }

	modifiedBogackiShampine23(const modifiedBogackiShampine23&);
   modifiedBogackiShampine23& operator=(const modifiedBogackiShampine23&);
   
   void modifiedRightHandSide(const G4double yIn[], G4double dydx[]);


 private:

	   G4double *ak2, *ak3, *ak4, *yTemp, *yIn;
      // for storing intermediate 'k' values in stepper
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations

    modifiedBogackiShampine23* fAuxStepper;


};

#endif /* modifiedBogackiShampine23 */





























