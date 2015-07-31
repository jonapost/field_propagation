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


#include "MagEqRhsbyTimestoreB.hh"


#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"

#include "G4ThreeVector.hh"

#include <iostream>
using namespace std;

MagEqRhs_byTime_storeB::MagEqRhs_byTime_storeB( G4MagneticField* MagField )
   : G4Mag_EqRhs( MagField ) {

   lastB_eval = new G4double[3];
}

MagEqRhs_byTime_storeB::~MagEqRhs_byTime_storeB() {
   delete[] lastB_eval;
}

void
MagEqRhs_byTime_storeB::EvaluateRhsGivenB( const G4double y[],
                      const G4double B[3],
                       G4double dydx[] ) const
{

   for (int i = 0; i < 3; i ++){
      lastB_eval[i] = B[i];
   }

   G4double cof = FCof();
   cof /= FMass();   // Correction because we are using dynamic variables
                     // position and velocity (not momentum).

   dydx[0] = y[3];
   dydx[1] = y[4];
   dydx[2] = y[5];

   dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)

   return ;
}

void
MagEqRhs_byTime_storeB::
 SetChargeMomentumMass( G4ChargeState particleCharge,
                        G4double MomentumXc,
         G4double mass)

{
   G4Mag_EqRhs::SetChargeMomentumMass( particleCharge, MomentumXc, mass);
}
