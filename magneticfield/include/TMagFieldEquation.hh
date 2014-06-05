


#include "G4ChargeState.hh"
#include "G4Mag_UsualEqRhs.hh"


template 
<class T_Field, class T_State>
class TMagFieldEquation : public G4Mag_UsualEqRhs
{
	public:

		TMagFieldEquation(T_Field* f)
			:	G4Mag_UsualEqRhs(f)
		{
			itsField = f;
		}

		~TMagFieldEquation(){;}

		inline void TGetFieldValue(const G4double Point[4],
				G4double Field[]) const
		{
			itsField->T_Field::GetFieldValue(Point, Field);
		}
		inline void TRightHandSide(T_State s) const
		{
			G4double Field[G4maximum_number_of_field_components]; 
			G4double  PositionAndTime[4];
			PositionAndTime[0] = s.x[0];
			PositionAndTime[1] = s.x[1];
			PositionAndTime[2] = s.x[2];
			PositionAndTime[3] = s.x[7];     
			GetFieldValue(PositionAndTime, Field) ;
			TEvaluateRhsGivenB(s, Field);
		}
		
		void TEvaluateRhsGivenB(T_State s, 
				const G4double B[3]) const;

	private:
		enum { G4maximum_number_of_field_components = 24 };
		T_Field* itsField;
};

