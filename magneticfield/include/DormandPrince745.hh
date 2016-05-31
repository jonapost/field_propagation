//  An implementation of the 5th order embedded RK method from the paper
//    J. R. Dormand and P. J. Prince, “A family of embedded Runge-Kutta formulae,”
//	    Journal of computational and applied …, vol. 6, no. 1, pp. 19–26, 1980.
//
//  DormandPrince7 - 5(4) embedded RK method
//
//  Design & Implementation by Somnath Banerjee
//  Supervision & code review: John Apostolakis
//
// Work supported by the Google Summer of Code 2015.
//
// First version: 25 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//
//  History
// ------------------------------------------
//  Created   : 25 May 2015.             - Somnath
//  Revisions : 
//   * 29 June 2015:  Added interpolate() method(s) - Somnath

#ifndef Dormand_Prince_745
#define Dormand_Prince_745

#include "G4MagIntegratorStepper.hh"

class DormandPrince745 : public G4MagIntegratorStepper
{
  public:
    DormandPrince745(G4EquationOfMotion *EqRhs,
					 G4int numberOfVariables = 6,
					 G4bool primary =  true);
    ~DormandPrince745();
   
    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    //For Preparing the Interpolant and calculating the extra stages
    void SetupInterpolation_low( /* const G4double yInput[],
                                    const G4double dydx[],
                                    const G4double Step */ );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate_low( /* const G4double yInput[],
                             const G4double dydx[],
                             const G4double Step, */ 
                          G4double yOut[],
                          G4double tau );
    
    inline void SetupInterpolation()
                            /* ( const G4double yInput[],
                                 const G4double dydx[],
                                 const G4double Step ) */
    { 
       SetupInterpolation_low( /* yInput, dydx, Step */ );
       // SetupInterpolation_high( /* yInput, dydx, Step */ );       
    }
    
    //For calculating the output at the tau fraction of Step
    inline void Interpolate(
                         /* const G4double yInput[],
                            const G4double dydx[],
                            const G4double Step,  */
                                  G4double tau,    
                                  G4double yOut[]
       )
    {
        Interpolate_low(  /* yInput, dydx, Step, */  yOut, tau);
        // Interpolate_high(  /* yInput, dydx, Step, */  yOut, tau);        
    }
    
    void SetupInterpolation_high( /* const G4double yInput[],
                               const G4double dydx[],
                               const G4double Step */ );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate_high( /* const G4double yInput[],
                              const G4double dydx[],
                              const G4double Step, */ 
                                 G4double yOut[],
                                 G4double tau );

    G4double  DistChord() const;
    G4double DistChord2() const;
    G4double DistChord3() const;
   
    //  Enabling method, with common code between implementations (and steppers)
    G4double DistLine( G4double yStart[], G4double yMid[], G4double yEnd[] ) const;   
    G4int IntegratorOrder() const {return 4; }
    
    //New copy constructor
    //  DormandPrince745(const DormandPrince745 &);
    
private :
    
    DormandPrince745& operator=(const DormandPrince745&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7,
      *ak8, *ak9, 	//For additional stages in the interpolant
      *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fInitialDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    DormandPrince745* fAuxStepper;
};
#endif /* defined(__Geant4__DormandPrince745__) */
