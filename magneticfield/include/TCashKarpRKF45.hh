//template mode of RKF45

#include "G4CashKarpRKF45.hh"

template
<class T_Equation, class T_Field, class T_State>
class TCashKarpRKF45 : public G4CashKarpRKF45
{
	public:
		
		TCashKarpRKF45(T_Equation* EqRhs)
			: G4CashKarpRKF45(EqRhs)
		{;}

		~TCashKarpRKF45(){;}

		void TStepper(T_State s,
				G4double h,
				G4double yerr[]);
};
