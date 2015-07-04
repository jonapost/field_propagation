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

class FineRKNG34: public G4MagIntegratorStepper {
public:
   virtual ~FineRKNG34();
   FineRKNG34(G4EquationOfMotion *EqRhs,
         G4int numberOfVariables = 6,
         G4bool primary = true);

   void Stepper( const G4double y[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[] ) ;

   G4double  DistChord()   const;
   G4int IntegratorOrder() const {return 4; }

   void InterpolatePosition(G4double xi, G4double yout[]);

private:

   G4double **a, **aprime;
   G4double **f;

   G4double b[5], bprime[5], c[5];
   G4double b_error[5], bprime_error[5];

   Interpolant *position_interpolant;

   G4double yInitial[8], yNext[3];
   G4double fInitial[8], fNext[3];

   //G4double last_step_len;
   // G4double last_time_value; // Hack to implement FSAL

   FineRKNG34 *fAuxStepper;

   G4double fLastStepLength;
   G4double *fLastInitialVector, *fLastFinalVector,
            *fLastDyDx, *fMidVector, *fMidError;
   // for DistChord calculations



};

#endif /* MAGNETICFIELD_INCLUDE_FINERKNG34_HH_ */
