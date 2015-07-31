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


   G4double c1 = 0./1., c2 = 1./10., c3 = 1./10., c4 = 3./10., c5 = 7./10., c6 = 7./10., c7 = 17./25., c8 = 17./25., c9 = 1./1.;

   G4double
   a21 = 1./10.,
   a31 = 1./20., a32 = 1./20.,
   a41 = -3./20., a42 = -27./20., a43 = 9./5.,
   a51 = 596./495., a52 = -4./11., a53 = -427./264., a54 = 5849./3960.,
   a61 = -30183643./75879540., a62 = -2893317./88748., a63 = 50178653./1264659., a64 = -40900531./18969885., a65 = 5./19.,
   a71 = -396583./95400., a72 = -169371./5000., a73 = 21382771./530000., a74 = -1735451./954000., a75 = 1./6., a76 = 1./53.,
   a81 = -7118348801./1589316300., a82 = -169371./5000., a83 = 4295617784819./105954420000., a84 = -160583069351./127145304000., a85 = -1292153./2665520., a86 = 3399687081./565090240., a87 = -17500./3029.,
   a91 = 2243113./2686068., a92 = 245./24., a93 = -19392535./1745568., a94 = 165115./361152., a95 = 193763./351120., a96 = 34289./24640., a97 = -461732./327845., a98 = 24232./327845.;

   G4double
   alpha21 = 1./200.,
   alpha31 = 1./200., alpha32 = 0./1.,
   alpha41 = -1./2200., alpha42 = 0./1., alpha43 = 1./22.,
   alpha51 = 637./6600., alpha52 = 0./1., alpha53 = -7./110., alpha54 = 7./33.,
   alpha61 = 637./6600., alpha62 = 0./1., alpha63 = -7./110., alpha64 = 7./33., alpha65 = 0./1.,
   alpha71 = 225437./1968750., alpha72 = 0./1., alpha73 = -30073./281250., alpha74 = 65569./281250., alpha75 = -5534276509./402334453125., alpha76 = 11940527048./2816341171875.,
   alpha81 = 225437./1968750., alpha82 = 0./1., alpha83 = -30073./281250., alpha84 = 65569./281250., alpha85 = -5534276509./402334453125., alpha86 = 11940527048./2816341171875., alpha87 = 0./1.,
   alpha91 = 151./2142., alpha92 = 0./1., alpha93 = 5./116., alpha94 = 385./1368., alpha95 = 6./55., alpha96 = 2017./9240., alpha97 = -8./35., alpha98 = 6058./983535.;

   G4double
   b1 = 151./2142., b2 = 0./1., b3 = 25./522., b4 = 275./684., b5 = 4./11., b6 = 2017./2772., b7 = -5./7., b8 = 15145./786828., b9 = 1./12.;

   G4double
   bhat1 = 1349./157500., bhat2 = 0./1., bhat3 = 7873./45000., bhat4 = 27457./90000., bhat5 = 130803094867./3218675., bhat6 = 1266768679133./62500., bhat7 = 0./3004097250000., bhat8 = -2./1., bhat9 = 11./5., bhat10 = 1./150.;

};

#endif /* MAGNETICFIELD_INCLUDE_MURUARKN6459_HH_ */
