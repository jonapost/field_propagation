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
#include "MuruaRKN6459.hh"

#include "G4Mag_EqRhs.hh"
#include "MagEqRhsbyTimestoreB.hh"


MuruaRKN6459::MuruaRKN6459(G4EquationOfMotion *EqRhs,
      G4int numberOfVariables): G4MagIntegratorStepper(EqRhs, numberOfVariables) {
   MagEqRhs_byTime_storeB *myField_as_storeBfield = dynamic_cast<MagEqRhs_byTime_storeB*>(GetEquationOfMotion());
}

MuruaRKN6459::~MuruaRKN6459() {
   // Don't delete B, MuruaRKN6459 doesn't own it.
}


void MuruaRKN6459::ComputeRhsWithStoredB( const G4double y[],
                                          //const G4double B[3],
                                          G4double dydx[] )
{
   G4double cof = myField_as_storeBfield -> FCof();
   cof /= myField_as_storeBfield -> FMass();   // Correction because we are using dynamic variables
                                // position and velocity (not momentum).
                                // (This is just a copy from Mag_UsualEqRhs_IntegrateByTime.)

   B = myField_as_storeBfield -> get_last_B_eval();

   // Ignore 0 through 3 components. We don't use them for a Nystrom stepper.
   /*
   dydx[0] = y[3];
   dydx[1] = y[4];
   dydx[2] = y[5];
   */

   dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)

   return ;
}


G4double  MuruaRKN6459::DistChord() const { // Dummy implementation

   return -1.;

}



void MuruaRKN6459::Stepper(  const G4double y[],
                  const G4double dydx[],
                  G4double h,
                  G4double yout[],
                  G4double yerr[] ) {

   int k;

   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c2*v0[k] + h*h*(alpha21*Vdot1[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a21*Vdot1[k+3]);
   }
   ComputeRightHandSide(ytemp, Vdot2);
   // Duplicate!: using last Q_i value
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a31*Vdot1[k+3] + a32*Vdot2[k+3]);
   }
   // Use last B field Evaluation
   B = myField_as_storeBfield -> get_last_B_eval();
   ComputeRhsWithStoredB(B, Vdot3);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c4*v0[k] + h*h*(alpha41*Vdot1[k+3]+alpha42*Vdot2[k+3]+alpha43*Vdot3[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a41*Vdot1[k+3] + a42*Vdot2[k+3] + a43*Vdot3[k+3]);
   }
   ComputeRightHandSide(ytemp, Vdot4);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c5*v0[k] + h*h*(alpha51*Vdot1[k+3]+alpha52*Vdot2[k+3]+alpha53*Vdot3[k+3]+alpha54*Vdot4[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a51*Vdot1[k+3] + a52*Vdot2[k+3] + a53*Vdot3[k+3] + a54*Vdot4[k+3]);
   }
   ComputeRightHandSide(ytemp, Vdot5);
   // Duplicate!: using last Q_i value
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a61*Vdot1[k+3] + a62*Vdot2[k+3] + a63*Vdot3[k+3] + a64*Vdot4[k+3] + a65*Vdot5[k+3]);
   }
   // Use last B field Evaluation
   B = myField_as_storeBfield -> get_last_B_eval();
   ComputeRhsWithStoredB(B, Vdot6);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c7*v0[k] + h*h*(alpha71*Vdot1[k+3]+alpha72*Vdot2[k+3]+alpha73*Vdot3[k+3]+alpha74*Vdot4[k+3]+alpha75*Vdot5[k+3]+alpha76*Vdot6[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a71*Vdot1[k+3] + a72*Vdot2[k+3] + a73*Vdot3[k+3] + a74*Vdot4[k+3] + a75*Vdot5[k+3] + a76*Vdot6[k+3]);
   }
   ComputeRightHandSide(ytemp, Vdot7);
   // Duplicate!: using last Q_i value
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a81*Vdot1[k+3] + a82*Vdot2[k+3] + a83*Vdot3[k+3] + a84*Vdot4[k+3] + a85*Vdot5[k+3] + a86*Vdot6[k+3] + a87*Vdot7[k+3]);
   }
   // Use last B field Evaluation
   B = myField_as_storeBfield -> get_last_B_eval();
   ComputeRhsWithStoredB(B, Vdot8);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c9*v0[k] + h*h*(alpha91*Vdot1[k+3]+alpha92*Vdot2[k+3]+alpha93*Vdot3[k+3]+alpha94*Vdot4[k+3]+alpha95*Vdot5[k+3]+alpha96*Vdot6[k+3]+alpha97*Vdot7[k+3]+alpha98*Vdot8[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a91*Vdot1[k+3] + a92*Vdot2[k+3] + a93*Vdot3[k+3] + a94*Vdot4[k+3] + a95*Vdot5[k+3] + a96*Vdot6[k+3] + a97*Vdot7[k+3] + a98*Vdot8[k+3]);
   }
   ComputeRightHandSide(ytemp, Vdot9);



}
