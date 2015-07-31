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
#ifndef MAGNETICFIELD_INCLUDE_MAGEQRHSBYTIMESTOREB_HH_
#define MAGNETICFIELD_INCLUDE_MAGEQRHSBYTIMESTOREB_HH_

#include "G4Mag_EqRhs.hh"



#include "G4Mag_EqRhs.hh"
#include "G4ChargeState.hh"

class MagEqRhs_byTime_storeB: public G4Mag_EqRhs {
public:
   MagEqRhs_byTime_storeB(G4MagneticField* MagField );

   virtual ~MagEqRhs_byTime_storeB();

   void EvaluateRhsGivenB(   const G4double y[],
                             const G4double B[3],
                                   G4double dydx[] ) const;
       // Given the value of the magnetic field B, this function
       // calculates the value of the derivative dydx.

   virtual void SetChargeMomentumMass( G4ChargeState particleCharge,
                                            G4double MomentumXc,
                                            G4double mass);


   inline G4double * get_last_B_eval() { return lastB_eval; }

//private:
   G4double *lastB_eval;

};



#endif /* MAGNETICFIELD_INCLUDE_MAGEQRHSBYTIMESTOREB_HH_ */
