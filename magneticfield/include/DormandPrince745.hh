//
//  DormandPrince745.hh
//  Geant4
//
//  Created by Somnath on 25/05/15.
//
//

#ifndef Dormand_Prince_745
#define Dormand_Prince_745

#include "G4MagIntegratorStepper.hh"

class DormandPrince745 : public G4MagIntegratorStepper
{
public:
    //Constructor using Equation
	DormandPrince745(G4EquationOfMotion *EqRhs,
					 G4int numberOfVariables = 6,
					 G4bool primary =  true);
	~DormandPrince745();

    //The Stepper
	void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 4; }
    
private :
    
    DormandPrince745(const DormandPrince745&);
        //Copy constructor kept private
    DormandPrince745& operator=(const DormandPrince745&); 
        //assignment operator overloaded
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    DormandPrince745* fAuxStepper;
};







#endif /* defined(__Geant4__DormandPrince745__) */
