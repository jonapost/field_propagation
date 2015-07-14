/*
 * FineRKNG45.hh
 *
 *  Created on: Jul 13, 2015
 *      Author: jason
 */

#ifndef MAGNETICFIELD_INCLUDE_FINERKNG45_HH_
#define MAGNETICFIELD_INCLUDE_FINERKNG45_HH_

#include "G4MagIntegratorStepper.hh"

class FineRKNG45: public G4MagIntegratorStepper {
public:
   virtual ~FineRKNG45();
   FineRKNG45(G4EquationOfMotion *EqRhs,
         G4int numberOfVariables = 6,
         G4bool primary = true);

   virtual void Stepper( const G4double y[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[] ) ;

   G4double  DistChord()   const;
   G4int IntegratorOrder() const {return 5; }

   //void InterpolatePosition(G4double xi, G4double yout[]);

private:

   G4double **a, **aprime;
   G4double **f;

   G4double b[7], bprime[7], c[7];
   G4double b_error[7], bprime_error[7];

   //Interpolant *position_interpolant;

   G4double yInitial[8], yNext[3];
   G4double fInitial[8], fNext[3];

   //G4double last_step_len;
   // G4double last_time_value; // Hack to implement FSAL

   FineRKNG45 *fAuxStepper;

   G4double fLastStepLength;
   G4double *fLastInitialVector, *fLastFinalVector,
            *fLastDyDx, *fMidVector, *fMidError;
   // for DistChord calculations



};

#endif /* MAGNETICFIELD_INCLUDE_FINERKNG45_HH_ */
