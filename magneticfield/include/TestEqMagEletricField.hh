#ifndef TESTEQMAGELECTRICFIELD_hh
#define TESTEQMAGELECTRICFIELD_hh

#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"
#include "G4ElectroMagneticField.hh"

template
<class State, class Field>
class TestEqMagEletricField : public G4EquationOfMotion
{
public:
	TestEqMagEletricField(Field *emField)
		: G4EquationOfMotion( emField ) {;}

	~TestEqMagEletricField() {;}

	void  SetChargeMomentumMass(G4ChargeState particleCharge,  
			G4double MomentumXc,
			G4double mass);

	//void EvaluateRhsGivenB(State s, const G4double Field[]); 
	void EvaluateRhsGivenB(const G4double y[], 
			G4double dydx[], const G4double Field[]) const;
        void TestEvaluateRhsGivenB(State s, G4double Field[]) const;

private:

	G4double        fElectroMagCof;
	G4double        fMassCof;
};

#endif

