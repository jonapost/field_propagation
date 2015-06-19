/*
 * ChawlaSharmaWrapper.hh
 *
 *  Created on: Jun 17, 2015
 *      Author: jason
 */

#ifndef MAGNETICFIELD_INCLUDE_CHAWLASHARMAWRAPPER_HH_
#define MAGNETICFIELD_INCLUDE_CHAWLASHARMAWRAPPER_HH_

#include "ChawlaSharmaRKNstepper.hh"

class ChawlaSharmaWrapper: public ChawlaSharmaRKNstepper {

public:
	ChawlaSharmaWrapper(G4Mag_EqRhs *EquationRhs,
            			G4int numberOfVariables = 6,
						G4int numStateVariables = 12);
	virtual ~ChawlaSharmaWrapper();

	void Stepper(const G4double yInput[],
				 const G4double dydx[],
				 	   G4double hstep,
					   G4double yOutput[],
					   G4double yError [] );
	void DumbStepper( const G4double y[],
						//const G4double dydx[],
						G4double h,
						G4double yOut[] );


	//void ComputeRightHandSide(const G4double y[], G4double dydx[] );

};

inline
ChawlaSharmaWrapper::ChawlaSharmaWrapper(G4Mag_EqRhs *EquationRhs,
                                     	 G4int numberOfVariables,
									 	 G4int numStateVariables)
     : ChawlaSharmaRKNstepper(EquationRhs, numberOfVariables, numStateVariables)
  { }





#endif /* MAGNETICFIELD_INCLUDE_CHAWLASHARMAWRAPPER_HH_ */
