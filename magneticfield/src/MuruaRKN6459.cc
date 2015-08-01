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
   myField_as_storeBfield = dynamic_cast<MagEqRhs_byTime_storeB*>(GetEquationOfMotion());



   c2 = 1./10., c3 = 1./10., c4 = 3./10., c5 = 7./10., c6 = 7./10., c7 = 17./25., c8 = 17./25., c9 = 1./1.;

   a21 = 1./10.,
   a31 = 1./20., a32 = 1./20.,
   a41 = -3./20., a42 = -27./20., a43 = 9./5.,
   a51 = 596./495., a52 = -4./11., a53 = -427./264., a54 = 5849./3960.,
   a61 = -30183643./75879540., a62 = -2893317./88748., a63 = 50178653./1264659., a64 = -40900531./18969885., a65 = 5./19.,
   a71 = -396583./95400., a72 = -169371./5000., a73 = 21382771./530000., a74 = -1735451./954000., a75 = 1./6., a76 = 1./53.,
   a81 = -7118348801./1589316300., a82 = -169371./5000., a83 = 4295617784819./105954420000., a84 = -160583069351./127145304000., a85 = -1292153./2665520., a86 = 3399687081./565090240., a87 = -17500./3029.,
   a91 = 2243113./2686068., a92 = 245./24., a93 = -19392535./1745568., a94 = 165115./361152., a95 = 193763./351120., a96 = 34289./24640., a97 = -461732./327845., a98 = 24232./327845.;



   alpha21 = 1./200.,
   alpha31 = 1./200., alpha32 = 0./1.,
   alpha41 = -1./2200., alpha42 = 0./1., alpha43 = 1./22.,
   alpha51 = 637./6600., alpha52 = 0./1., alpha53 = -7./110., alpha54 = 7./33.,
   alpha61 = 637./6600., alpha62 = 0./1., alpha63 = -7./110., alpha64 = 7./33., alpha65 = 0./1.,
   alpha71 = 225437./1968750., alpha72 = 0./1., alpha73 = -30073./281250., alpha74 = 65569./281250., alpha75 = -5534276509./402334453125., alpha76 = 11940527048./2816341171875.,
   alpha81 = 225437./1968750., alpha82 = 0./1., alpha83 = -30073./281250., alpha84 = 65569./281250., alpha85 = -5534276509./402334453125., alpha86 = 11940527048./2816341171875., alpha87 = 0./1.,
   alpha91 = 151./2142., alpha92 = 0./1., alpha93 = 5./116., alpha94 = 385./1368., alpha95 = 6./55., alpha96 = 2017./9240., alpha97 = -8./35., alpha98 = 6058./983535.;

   b1 = 151./2142., b2 = 0./1., b3 = 25./522., b4 = 275./684., b5 = 4./11., b6 = 2017./2772., b7 = -5./7., b8 = 15145./786828., b9 = 1./12.;
   bhat1 = 1349./157500., bhat2 = 0./1., bhat3 = 7873./45000., bhat4 = 27457./90000., bhat5 = 130803094867./3218675., bhat6 = 1266768679133./62500., bhat7 = 0./3004097250000., bhat8 = -2./1., bhat9 = 11./5., bhat10 = 1./150.;
   beta1 = 151./2142., beta2 = 0./1., beta3 = 5./116., beta4 = 385./1368., beta5 = 6./55., beta6 = 2017./9240., beta7 = -8./35., beta8 = 6058./983535., beta9 = 0./1.;
   beta_hat1 = 1349./157500., beta_hat2 = 0./1., beta_hat3 = 7873./50000., beta_hat4 = 192199./900000., beta_hat5 = 130803094867./1072891875000., beta_hat6 = 1266768679133./10013657500000., beta_hat7 = 0./1., beta_hat8 = -16./125., beta_hat9 = 0./1., beta_hat10 = 0./1.;


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
      q0[k] = y[k];
      v0[k] = y[k + 3];
      Vdot1[k + 3] = dydx[k + 3];
   }


   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c2*v0[k] + h*h*(alpha21*Vdot1[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a21*Vdot1[k+3]);
   }
   RightHandSide(ytemp, Vdot2);
   // Duplicate!: using last Q_i value
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a31*Vdot1[k+3] + a32*Vdot2[k+3]);
   }
   // Use last B field Evaluation
   B = myField_as_storeBfield -> get_last_B_eval();
   ComputeRhsWithStoredB(ytemp, Vdot3);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c4*v0[k] + h*h*(alpha41*Vdot1[k+3]+alpha42*Vdot2[k+3]+alpha43*Vdot3[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a41*Vdot1[k+3] + a42*Vdot2[k+3] + a43*Vdot3[k+3]);
   }
   RightHandSide(ytemp, Vdot4);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c5*v0[k] + h*h*(alpha51*Vdot1[k+3]+alpha52*Vdot2[k+3]+alpha53*Vdot3[k+3]+alpha54*Vdot4[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a51*Vdot1[k+3] + a52*Vdot2[k+3] + a53*Vdot3[k+3] + a54*Vdot4[k+3]);
   }
   RightHandSide(ytemp, Vdot5);
   // Duplicate!: using last Q_i value
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a61*Vdot1[k+3] + a62*Vdot2[k+3] + a63*Vdot3[k+3] + a64*Vdot4[k+3] + a65*Vdot5[k+3]);
   }
   // Use last B field Evaluation
   B = myField_as_storeBfield -> get_last_B_eval();
   ComputeRhsWithStoredB(ytemp, Vdot6);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c7*v0[k] + h*h*(alpha71*Vdot1[k+3]+alpha72*Vdot2[k+3]+alpha73*Vdot3[k+3]+alpha74*Vdot4[k+3]+alpha75*Vdot5[k+3]+alpha76*Vdot6[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a71*Vdot1[k+3] + a72*Vdot2[k+3] + a73*Vdot3[k+3] + a74*Vdot4[k+3] + a75*Vdot5[k+3] + a76*Vdot6[k+3]);
   }
   RightHandSide(ytemp, Vdot7);
   // Duplicate!: using last Q_i value
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a81*Vdot1[k+3] + a82*Vdot2[k+3] + a83*Vdot3[k+3] + a84*Vdot4[k+3] + a85*Vdot5[k+3] + a86*Vdot6[k+3] + a87*Vdot7[k+3]);
   }
   // Use last B field Evaluation
   B = myField_as_storeBfield -> get_last_B_eval();
   ComputeRhsWithStoredB(ytemp, Vdot8);
   for (k = 0; k < 3; k ++) {
      ytemp[k] = q0[k] + h*c9*v0[k] + h*h*(alpha91*Vdot1[k+3]+alpha92*Vdot2[k+3]+alpha93*Vdot3[k+3]+alpha94*Vdot4[k+3]+alpha95*Vdot5[k+3]+alpha96*Vdot6[k+3]+alpha97*Vdot7[k+3]+alpha98*Vdot8[k+3]);
   }
   for (k = 0; k < 3; k ++) {
      ytemp[k + 3] = v0[k]+h*(a91*Vdot1[k+3] + a92*Vdot2[k+3] + a93*Vdot3[k+3] + a94*Vdot4[k+3] + a95*Vdot5[k+3] + a96*Vdot6[k+3] + a97*Vdot7[k+3] + a98*Vdot8[k+3]);
   }
   RightHandSide(ytemp, Vdot9);

   G4double temp_var;
   for (k = 3; k < 6; k ++) {
      temp_var = beta1*Vdot1[k] + beta2*Vdot2[k] + beta3*Vdot3[k] + beta4*Vdot4[k] + beta5*Vdot5[k] + beta6*Vdot6[k] + beta7*Vdot7[k] + beta8*Vdot8[k] + beta9*Vdot9[k];
      yout[k - 3] = q0[k - 3] + h*v0[k - 3] + h*h*temp_var;
   }
   for (k = 3; k < 6; k ++) {
      temp_var = b1*Vdot1[k] + b2*Vdot2[k] + b3*Vdot3[k] + b4*Vdot4[k] + b5*Vdot5[k] + b6*Vdot6[k] + b7*Vdot7[k] + b8*Vdot8[k] + b9*Vdot9[k];
      yout[k] = v0[k - 3] + h*temp_var;
   }


}
