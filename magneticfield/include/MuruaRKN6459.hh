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
#ifndef MAGNETICFIELD_INCLUDE_MURUARKN6459_HH_
#define MAGNETICFIELD_INCLUDE_MURUARKN6459_HH_


#define NO_STATE_VARIABLES 12
#define NO_INTEGRATION_VARIABLES 6



#include "G4MagIntegratorStepper.hh"
#include "MagEqRhsbyTimestoreB.hh"

class MuruaRKN6459: public G4MagIntegratorStepper {
public:
   MuruaRKN6459(G4EquationOfMotion *EqRhs,
         G4int numberOfVariables);
   virtual ~MuruaRKN6459();

   virtual void Stepper(  const G4double y[],
                  const G4double dydx[],
                  G4double h,
                  G4double yout[],
                  G4double yerr[] ) ;

   void ComputeRhsWithStoredB(const G4double y[],
         //const G4double B[3],
         G4double Rhs[] );

   G4double  DistChord()   const;
   G4int IntegratorOrder() const {return 6; }

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
   G4double
   bhat1, bhat2, bhat3, bhat4, bhat5, bhat6, bhat7, bhat8, bhat9, bhat10;
   G4double
   beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9;
   G4double
   beta_hat1, beta_hat2, beta_hat3, beta_hat4, beta_hat5, beta_hat6, beta_hat7, beta_hat8, beta_hat9, beta_hat10;

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
            Vdot9[NO_INTEGRATION_VARIABLES];









   /*


   G4double
   c2 = 1./10., c3 = 1./10., c4 = 3./10., c5 = 7./10., c6 = 7./10., c7 = 17./25., c8 = 17./25., c9 = 1./1.;

   G4double
   a21 = 1./10.,
   a31 = 1./20., a32 = 1./20.,
   a41 = -3./20., a42 = -27./20., a43 = 9./5.,
   a51 = 596./495., a52 = -4./11., a53 = -427./264., a54 = 5849./3960.,
   a61 = -30183643./75879540., a62 = -2893317./88748., a63 = 50178653./1264659., a64 = -40900531./18969885., a65 = 5./19.,
   a71 = -396583./95400., a72 = -169371./5000., a73 = 21382771./530000., a74 = -1735451./954000., a75 = 1./6., a76 = 1./53.,
   a81 = -7118348801./1589316300., a82 = -169371./5000., a83 = 4295617784819./105954420000., a84 = -160583069351./127145304000., a85 = -1292153./2665520., a86 = 3399687081./565090240., a87 = -17500./3029.,
   a91 = 2243113./2686068., a92 = 245./24., a93 = -19392535./1745568., a94 = 165115./361152., a95 = 193763./351120., a96 = 34289./24640., a97 = -461732./327845., a98 = 24232./327845.;

   const static G4double
   alpha21 = 1./200.,
   alpha31 = 1./200., alpha32 = 0./1.,
   alpha41 = -1./2200., alpha42 = 0./1., alpha43 = 1./22.,
   alpha51 = 637./6600., alpha52 = 0./1., alpha53 = -7./110., alpha54 = 7./33.,
   alpha61 = 637./6600., alpha62 = 0./1., alpha63 = -7./110., alpha64 = 7./33., alpha65 = 0./1.,
   alpha71 = 225437./1968750., alpha72 = 0./1., alpha73 = -30073./281250., alpha74 = 65569./281250., alpha75 = -5534276509./402334453125., alpha76 = 11940527048./2816341171875.,
   alpha81 = 225437./1968750., alpha82 = 0./1., alpha83 = -30073./281250., alpha84 = 65569./281250., alpha85 = -5534276509./402334453125., alpha86 = 11940527048./2816341171875., alpha87 = 0./1.,
   alpha91 = 151./2142., alpha92 = 0./1., alpha93 = 5./116., alpha94 = 385./1368., alpha95 = 6./55., alpha96 = 2017./9240., alpha97 = -8./35., alpha98 = 6058./983535.;

   const static G4double
   b1 = 151./2142., b2 = 0./1., b3 = 25./522., b4 = 275./684., b5 = 4./11., b6 = 2017./2772., b7 = -5./7., b8 = 15145./786828., b9 = 1./12.;

   const static G4double
   bhat1 = 1349./157500., bhat2 = 0./1., bhat3 = 7873./45000., bhat4 = 27457./90000., bhat5 = 130803094867./3218675., bhat6 = 1266768679133./62500., bhat7 = 0./3004097250000., bhat8 = -2./1., bhat9 = 11./5., bhat10 = 1./150.;

   const static G4double
   beta1 = 151./2142., beta2 = 0./1., beta3 = 5./116., beta4 = 385./1368., beta5 = 6./55., beta6 = 2017./9240., beta7 = -8./35., beta8 = 6058./983535., beta9 = 0./1.;

   const static G4double
   beta_hat1 = 1349./157500., beta_hat2 = 0./1., beta_hat3 = 7873./50000., beta_hat4 = 192199./900000., beta_hat5 = 130803094867./1072891875000., beta_hat6 = 1266768679133./10013657500000., beta_hat7 = 0./1., beta_hat8 = -16./125., beta_hat9 = 0./1., beta_hat10 = 0./1.;

   */


};
/*
G4double
MuruaRKN6459::a21 = 1./10.,
MuruaRKN6459::a31 = 1./20., MuruaRKN6459::a32 = 1./20.,
MuruaRKN6459::a41 = -3./20., MuruaRKN6459::a42 = -27./20., MuruaRKN6459::a43 = 9./5.,
MuruaRKN6459::a51 = 596./495., MuruaRKN6459::a52 = -4./11., MuruaRKN6459::a53 = -427./264., MuruaRKN6459::a54 = 5849./3960.,
MuruaRKN6459::a61 = -30183643./75879540., MuruaRKN6459::a62 = -2893317./88748., MuruaRKN6459::a63 = 50178653./1264659., MuruaRKN6459::a64 = -40900531./18969885., MuruaRKN6459::a65 = 5./19.,
MuruaRKN6459::a71 = -396583./95400., MuruaRKN6459::a72 = -169371./5000., MuruaRKN6459::a73 = 21382771./530000., MuruaRKN6459::a74 = -1735451./954000., MuruaRKN6459::a75 = 1./6., MuruaRKN6459::a76 = 1./53.,
MuruaRKN6459::a81 = -7118348801./1589316300., MuruaRKN6459::a82 = -169371./5000., MuruaRKN6459::a83 = 4295617784819./105954420000., MuruaRKN6459::a84 = -160583069351./127145304000., MuruaRKN6459::a85 = -1292153./2665520., MuruaRKN6459::a86 = 3399687081./565090240., MuruaRKN6459::a87 = -17500./3029.,
MuruaRKN6459::a91 = 2243113./2686068., MuruaRKN6459::a92 = 245./24., MuruaRKN6459::a93 = -19392535./1745568., MuruaRKN6459::a94 = 165115./361152., MuruaRKN6459::a95 = 193763./351120., MuruaRKN6459::a96 = 34289./24640., MuruaRKN6459::a97 = -461732./327845., MuruaRKN6459::a98 = 24232./327845.;

G4double
MuruaRKN6459::alpha21,
MuruaRKN6459::alpha31, MuruaRKN6459::alpha32,
MuruaRKN6459::alpha41, MuruaRKN6459::alpha42, MuruaRKN6459::alpha43,
MuruaRKN6459::alpha51, MuruaRKN6459::alpha52, MuruaRKN6459::alpha53, MuruaRKN6459::alpha54,
MuruaRKN6459::alpha61, MuruaRKN6459::alpha62, MuruaRKN6459::alpha63, MuruaRKN6459::alpha64, MuruaRKN6459::alpha65,
MuruaRKN6459::alpha71, MuruaRKN6459::alpha72, MuruaRKN6459::alpha73, MuruaRKN6459::alpha74, MuruaRKN6459::alpha75, MuruaRKN6459::alpha76,
MuruaRKN6459::alpha81, MuruaRKN6459::alpha82, MuruaRKN6459::alpha83, MuruaRKN6459::alpha84, MuruaRKN6459::alpha85, MuruaRKN6459::alpha86, MuruaRKN6459::alpha87,
MuruaRKN6459::alpha91, MuruaRKN6459::alpha92, MuruaRKN6459::alpha93, MuruaRKN6459::alpha94, MuruaRKN6459::alpha95, MuruaRKN6459::alpha96, MuruaRKN6459::alpha97, MuruaRKN6459::alpha98;

G4double
MuruaRKN6459::b1 = 151./2142., MuruaRKN6459::b2 = 0./1., MuruaRKN6459::b3 = 25./522., MuruaRKN6459::b4 = 275./684., MuruaRKN6459::b5 = 4./11., MuruaRKN6459::b6 = 2017./2772., MuruaRKN6459::b7 = -5./7., MuruaRKN6459::b8 = 15145./786828., MuruaRKN6459::b9 = 1./12.;
G4double
MuruaRKN6459::bhat1 = 1349./157500., MuruaRKN6459::bhat2 = 0./1., MuruaRKN6459::bhat3 = 7873./45000., MuruaRKN6459::bhat4 = 27457./90000., MuruaRKN6459::bhat5 = 130803094867./3218675., MuruaRKN6459::bhat6 = 1266768679133./62500., MuruaRKN6459::bhat7 = 0./3004097250000., MuruaRKN6459::bhat8 = -2./1., MuruaRKN6459::bhat9 = 11./5., MuruaRKN6459::bhat10 = 1./150.;
G4double
MuruaRKN6459::beta1 = 151./2142., MuruaRKN6459::beta2 = 0./1., MuruaRKN6459::beta3 = 5./116., MuruaRKN6459::beta4 = 385./1368., MuruaRKN6459::beta5 = 6./55., MuruaRKN6459::beta6 = 2017./9240., MuruaRKN6459::beta7 = -8./35., MuruaRKN6459::beta8 = 6058./983535., MuruaRKN6459::beta9 = 0./1.;
G4double
MuruaRKN6459::beta_hat1 = 1349./157500., MuruaRKN6459::beta_hat2 = 0./1., MuruaRKN6459::beta_hat3 = 7873./50000., MuruaRKN6459::beta_hat4 = 192199./900000., MuruaRKN6459::beta_hat5 = 130803094867./1072891875000., MuruaRKN6459::beta_hat6 = 1266768679133./10013657500000., MuruaRKN6459::beta_hat7 = 0./1., MuruaRKN6459::beta_hat8 = -16./125., MuruaRKN6459::beta_hat9 = 0./1., MuruaRKN6459::beta_hat10 = 0./1.;

G4double
MuruaRKN6459::c2 = 1./10., MuruaRKN6459::c3 = 1./10., MuruaRKN6459::c4 = 3./10., MuruaRKN6459::c5 = 7./10., MuruaRKN6459::c6 = 7./10., MuruaRKN6459::c7 = 17./25., MuruaRKN6459::c8 = 17./25., MuruaRKN6459::c9 = 1./1.;


*/


#endif /* MAGNETICFIELD_INCLUDE_MURUARKN6459_HH_ */
