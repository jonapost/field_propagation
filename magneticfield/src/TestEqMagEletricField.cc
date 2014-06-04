#include "TestEqMagElectricField.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

void  
TestEqMagElectricField::SetChargeMomentumMass
(G4ChargeState particleCharge, G4double, G4double particleMass)
{
	G4double pcharge = particleCharge.GetCharge();
	fElectroMagCof =  eplus*pcharge*c_light ;
	fMassCof = particleMass*particleMass ; 
}

void 
TestEqMagElectricField::EvaluateRhsGivenB
(const G4double y[], G4double dydx[], const G4double Field[]) const
{
State s;
s.x = y;
TestEvaluateRhsGivenB(s, Field);
}


void
TestEqMagElectricField::TestEvaluateRhsGivenB
(State s, G4double Field[]) const 
{
	// Components of y:
	//    0-2 dr/ds, 
	//    3-5 dp/ds - momentum derivatives 

	G4double pSquared = s.x[3]*s.x[3] + s.x[4]*s.x[4] + s.x[5]*s.x[5] ;

	G4double Energy   = std::sqrt( pSquared + fMassCof );
	G4double cof2     = Energy/c_light ;

	G4double pModuleInverse  = 1.0/std::sqrt(pSquared) ;

	//  G4double inverse_velocity = Energy * c_light * pModuleInverse;
	G4double inverse_velocity = Energy * pModuleInverse / c_light;

	G4double cof1     = fElectroMagCof*pModuleInverse ;

	//  G4double vDotE = s.x[3]*Field[3] + s.x[4]*Field[4] + s.x[5]*Field[5] ;


	s.xp[0] = s.x[3]*pModuleInverse ;                         
	s.xp[1] = s.x[4]*pModuleInverse ;                         
	s.xp[2] = s.x[5]*pModuleInverse ;                        

	s.xp[3] = cof1*(cof2*Field[3] + (s.x[4]*Field[2] - s.x[5]*Field[1])) ;

	s.xp[4] = cof1*(cof2*Field[4] + (s.x[5]*Field[0] - s.x[3]*Field[2])) ; 

	s.xp[5] = cof1*(cof2*Field[5] + (s.x[3]*Field[1] - s.x[4]*Field[0])) ;  

	s.xp[6] = 0.;//not used

	// Lab Time of flight
	s.xp[7] = inverse_velocity;
	return ;
}
