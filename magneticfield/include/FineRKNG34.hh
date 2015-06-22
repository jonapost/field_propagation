/*
 * FineRKNG34.hh
 *
 *  Created on: Jun 21, 2015
 *      Author: jason
 */

#ifndef MAGNETICFIELD_INCLUDE_FINERKNG34_HH_
#define MAGNETICFIELD_INCLUDE_FINERKNG34_HH_

#include "G4MagIntegratorStepper.hh"

class FineRKNG34: public G4MagIntegratorStepper {
public:
   virtual ~FineRKNG34();
   FineRKNG34(G4EquationOfMotion *EqRhs,
         G4int numberOfVariables = 6,
         G4int numberOfStateVariables = 12);

   void Stepper( const G4double y[],
                    const G4double dydx[],
                    G4double h,
                    G4double yout[],
                    G4double yerr[] ) ;

   G4double  DistChord()   const;
   G4int IntegratorOrder() const {return 4; }

private:

   G4double **a, **aprime;
   G4double **f;

   G4double b[5], bprime[5], c[5];
   G4double b_error[5], bprime_error[5];




};

#endif /* MAGNETICFIELD_INCLUDE_FINERKNG34_HH_ */
