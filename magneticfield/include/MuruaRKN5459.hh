// Nystrom stepper implemenation by Jason Suagee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 27 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
#ifndef MAGNETICFIELD_INCLUDE_MURUARKN5459_HH_
#define MAGNETICFIELD_INCLUDE_MURUARKN5459_HH_

#define NO_STATE_VARIABLES 12
#define NO_INTEGRATION_VARIABLES 6



#include "G4MagIntegratorStepper.hh"
#include "MagEqRhsbyTimestoreB.hh"
#include "Interpolant.hh"


class MuruaRKN5459: public G4MagIntegratorStepper {
public:
   MuruaRKN5459(G4EquationOfMotion *EqRhs,
         G4int numberOfVariables,
         G4bool primary = true);
   virtual ~MuruaRKN5459();

   virtual void Stepper(  const G4double y[],
                  const G4double dydx[],
                  G4double h,
                  G4double yout[],
                  G4double yerr[] ) ;

   void ComputeRhsWithStoredB(const G4double y[],
         //const G4double B[3],
         G4double Rhs[] );

   G4double  DistChord()   const;
   G4int IntegratorOrder() const {return 5; }

   void set_MagEqRhs_storedBfield( G4EquationOfMotion *EqRhs );

private:

   G4double
   a21,
   a31, a32,
   a41, a42, a43,
   a51, a52, a53, a54,
   a61, a62, a63, a64, a65,
   a71, a72, a73, a74, a75, a76,
   a81, a82, a83, a84, a85, a86, a87,
   a91, a92, a93, a94, a95, a96, a97, a98;

   G4double
   alpha21,
   alpha31, alpha32,
   alpha41, alpha42, alpha43,
   alpha51, alpha52, alpha53, alpha54,
   alpha61, alpha62, alpha63, alpha64, alpha65,
   alpha71, alpha72, alpha73, alpha74, alpha75, alpha76,
   alpha81, alpha82, alpha83, alpha84, alpha85, alpha86, alpha87,
   alpha91, alpha92, alpha93, alpha94, alpha95, alpha96, alpha97, alpha98;

   G4double
   b1, b2, b3, b4, b5, b6, b7, b8, b9;

   //G4double
   //bhat1, bhat2, bhat3, bhat4, bhat5, bhat6, bhat7, bhat8, bhat9, bhat10;

   G4double
   berr1,  berr2,  berr3,  berr4,  berr5,  berr6,  berr7,  berr8,  berr9,  berr10;

   G4double
   beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9;

   //G4double
   //beta_hat1, beta_hat2, beta_hat3, beta_hat4, beta_hat5, beta_hat6, beta_hat7, beta_hat8, beta_hat9, beta_hat10;

   G4double
   beta_err1,  beta_err2,  beta_err3,  beta_err4,  beta_err5,  beta_err6,  beta_err7,  beta_err8,  beta_err9,  beta_err10;

   G4double
   c2, c3, c4, c5, c6, c7, c8, c9; // c1 is always 0 (and not even used) so we don't include it.




   MagEqRhs_byTime_storeB *myField_as_storeBfield;

   G4double *B;

   G4double q0[3], v0[3];
   G4double ytemp[NO_INTEGRATION_VARIABLES]; // Long enough just for ComputeRightHandSide()
   G4double Vdot1[NO_INTEGRATION_VARIABLES],
            Vdot2[NO_INTEGRATION_VARIABLES],
            Vdot3[NO_INTEGRATION_VARIABLES],
            Vdot4[NO_INTEGRATION_VARIABLES],
            Vdot5[NO_INTEGRATION_VARIABLES],
            Vdot6[NO_INTEGRATION_VARIABLES],
            Vdot7[NO_INTEGRATION_VARIABLES],
            Vdot8[NO_INTEGRATION_VARIABLES],
            Vdot9[NO_INTEGRATION_VARIABLES],
            Vdot10[NO_INTEGRATION_VARIABLES];


   Interpolant *position_interpolant;

   G4double yInitial[8], yNext[3];
   G4double fInitial[8], fNext[3];

   //G4double fLastDyDx[3], fNextDyDx[3];
   G4double fLastDyDx[6], fNextDyDx[6];      // Temporary (just changed to 6 to work with aux stepper).

   //G4double last_step_len;
   // G4double last_time_value; // Hack to implement FSAL

   //FineRKNG45 *fAuxStepper;

   G4double fLastStepLength;
   G4double *fLastInitialVector, *fLastFinalVector,
            *fMidVector, *fMidError;
   // for DistChord calculations



   MuruaRKN5459 *fAuxStepper;



};



#endif /* MAGNETICFIELD_INCLUDE_MURUARKN5459_HH_ */
