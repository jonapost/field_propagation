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

/*
 * Mag_UsualEqRhs_IntegrateByTime.hh
 *
 *  Created on: Jun 11, 2015
 *      Author: jason
 */

#ifndef MAGNETICFIELD_INCLUDE_MAG_USUALEQRHS_INTEGRATEBYTIME_HH_
#define MAGNETICFIELD_INCLUDE_MAG_USUALEQRHS_INTEGRATEBYTIME_HH_



#include "G4Mag_EqRhs.hh"
#include "G4ChargeState.hh"
// temp include:
#include <iomanip>


class Mag_UsualEqRhs_IntegrateByTime: public G4Mag_EqRhs {
public:
	Mag_UsualEqRhs_IntegrateByTime(G4MagneticField* MagField );

	virtual ~Mag_UsualEqRhs_IntegrateByTime();

	void EvaluateRhsGivenB( const G4double y[],
	                             const G4double B[3],
	                                   G4double dydx[] ) const;
	       // Given the value of the magnetic field B, this function
	       // calculates the value of the derivative dydx.

	virtual void SetChargeMomentumMass( G4ChargeState particleCharge,
	                                         G4double MomentumXc,
	                                         G4double mass);
};

#endif /* MAGNETICFIELD_INCLUDE_MAG_USUALEQRHS_INTEGRATEBYTIME_HH_ */
