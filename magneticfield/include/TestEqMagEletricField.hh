#ifndef TESTEQMAGELECTRICFIELD_hh
#define TESTEQMAGELECTRICFIELD_hh

#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"
#include "G4ElectroMagneticField.hh"

template
<class State, class Field>
class TestEqMagEletricField
{
public:
	TestEqMagEletricField(Field *emField)
		: G4EquationOfMotion( emField ) {;}

	~TestEqMagEletricField() {;}

	void  SetChargeMomentumMass(G4ChargeState particleCharge,  
			G4double MomentumXc,
			G4double mass);

	void EvaluateRhsGivenB(State s, const State f) const;

private:

	G4double        fElectroMagCof ;
	G4double        fMassCof;
};

#endif

