//experinmental purpose only


#ifndef G4_EquationOfMotion_DEF
#define G4_EquationOfMotion_DEF

#include "G4Types.hh"      // "globals.hh"
#include "G4Field.hh"
#include "G4ChargeState.hh"

template 
<class Equation>
class G4EquationOfMotion 
{
public:  

	G4EquationOfMotion(G4Field* field)
	{
		Equation itsEq(field);
	}

	virtual ~G4EquationOfMotion();

	inline void EvaluateRhsGivenB( const  G4double y[],
			const  G4double B[3],
			G4double dydx[] )
	{
		itsEq->EvaluateRhsGivenB(y, B, dydx);
	} 


	inline void SetChargeMomentumMass(G4ChargeState particleCharge,
			G4double MomentumXc,
			G4double MassXc2)
	{
		itsEq->SetChargeMomentumMass(particleCharge, 
				MomentumXc,
				MassXc2);
	} 

	inline void RightHandSide( const  G4double y[],
			G4double dydx[] ) const;

	void EvaluateRhsReturnB( const  G4double y[],
			G4double dydx[],
			G4double Field[]  ) const;

	inline void GetFieldValue( const  G4double Point[4],
			G4double Field[] )  const
	{
		itsEq->GetFieldValue(Point[4], Field[]);
	}

	/* to enable this part, we need to use: 
	 * template <class Equation, class Field> instead
	 inline     
	 const Field* 
	 GetFieldObj() const
	 {
	 return itsEq->GetFieldObj();
	 }
	 inline
	 void  SetFieldObj(Field* pField)
	 {
	 itsEq->SetFieldObj(pField);
	 }

	 */
private:
	// const int G4maximum_number_of_field_components = 24;
	enum { G4maximum_number_of_field_components = 24 } ;

	Equation *itsEq;

};

//#include "G4EquationOfMotion.icc"
#endif /* G4_EquationOfMotion_DEF */
