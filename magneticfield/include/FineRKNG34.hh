/*
 * FineRKNG34.hh
 *
 *  Created on: Jun 21, 2015
 *      Author: Jason Suagee
 */

#ifndef MAGNETICFIELD_INCLUDE_FINERKNG34_HH_
#define MAGNETICFIELD_INCLUDE_FINERKNG34_HH_

#include "G4MagIntegratorStepper.hh"
#include "Interpolant.hh"

#define NUMBER_INTERPOLATION_VARIABLES 6
#define NO_STATE_VARIABLES 12
#define NO_POSITION_VARIABLES 3



class FineRKNG34: public G4MagIntegratorStepper {
public:
   virtual ~FineRKNG34();
   FineRKNG34(G4EquationOfMotion *EqRhs,
         G4int numberOfVariables = 6,
         G4bool primary = true);

   virtual void Stepper( const G4double y[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[] ) ;

   G4double  DistChord()   const;
   G4int IntegratorOrder() const {return 3; }   // Was 4, but changed it to order of lower stepper
                                                // as per p. 168 of Hairer, Norset & Wanner (section on automatic step-size control).

   //void InterpolatePosition(G4double xi, G4double yout[]);

private:

   G4double **a, **aprime;
   G4double **f;

   G4double b[5], bprime[5], c[5];
   G4double b_error[5], bprime_error[5];

   Interpolant *position_interpolant;

   G4double yInitial[NO_STATE_VARIABLES], yNext[NO_STATE_VARIABLES];
   G4double yPrimeInitial[NO_STATE_VARIABLES], yPrimeNext[NO_STATE_VARIABLES];
   G4double fInitial[NO_POSITION_VARIABLES], fNext[NO_POSITION_VARIABLES];

   G4double fLastStepLength;

   // for DistChord calculations
   G4double *fMidVector;

};

#endif /* MAGNETICFIELD_INCLUDE_FINERKNG34_HH_ */
