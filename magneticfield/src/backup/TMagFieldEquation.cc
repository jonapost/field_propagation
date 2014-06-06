/*
#include "TMagFieldEquation.hh"

template 
<class T_Field>
void 
TMagFieldEquation<T_Field>::TEvaluateRhsGivenB(const G4double y[],
		const G4double B[3],
		G4double dydx[]) const
{
	G4double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
	G4double inv_momentum_magnitude = 1.0 / std::sqrt( momentum_mag_square );

	G4double cof = FCof()*inv_momentum_magnitude;

	dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
	dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
	dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

	dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
	dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
	dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)

	return ;
}
*/
