//  An implementation of the embedded RK method from the following paper 
//  by P. Bogacki and L. F. Shampine:
//    “An efficient Runge-Kutta (4,5) pair,” 
//   Comput. Math. with Appl., vol. 32, no. 6, pp. 15–28, Sep. 1996.
//
//  An interpolation method provides the value of an intermediate
//   point in a step -- if a step was sucessful.
//
//  This version can provide the FSAL property of the method,
//    which allows the reuse of the last derivative in the next step.
//    but only by using the additional method GetLastDyDx.
//  (Alternative interface for simpler use of FSAL is under development.)
//
//  Design & Implementation by Somnath Banerjee
//  Supervision & code review: John Apostolakis
//
// Work supported by the Google Summer of Code 2015.
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//
// Note: 
//   This version does not utilise the FSAL property of the method,
//   which would allow the reuse of the last derivative in the next step.
//   (Alternative FSAL implementation created with revised interface)
///////////////////////////////////////////////////////////////////////////////

#ifndef Bogacki_Shampine_45
#define Bogacki_Shampine_45

#include "G4MagIntegratorStepper.hh"

class BogackiShampine45 : public G4MagIntegratorStepper
{
public:
    BogackiShampine45(G4EquationOfMotion *EqRhs,
                        G4int numberOfVariables = 6,
                        G4bool primary =  true);
    ~BogackiShampine45();
    
    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    // This Stepper provides 'dense output'.  After a successful
    //  step, it is possible to obtain an estimate of the value
    //  of the function at an intermediate point of the interval.
    //  This requires only two additional evaluations of the
    //  derivative (and thus the field).

   inline void SetupInterpolation()
     //  const G4double yInput[], const G4double dydx[], const G4double Step )
    {
       SetupInterpolationHigh(); // ( yInput, dydx, Step);
    }
    
    //For calculating the output at the tau fraction of Step
    inline void Interpolate( // const G4double yInput[],
                             // const G4double dydx[],
                             // const G4double Step,
                             G4double tau, 
                             G4double yOut[]) // Output value
    {
       InterpolateHigh( tau, yOut);       
       // InterpolateHigh( yInput, dydx, Step, yOut, tau);
    }
    
    void SetupInterpolationHigh();
         // const G4double yInput[], const G4double, const G4double Step );
    
    // For calculating the output at the tau fraction of Step
    void InterpolateHigh( // const G4double yInput[],
                          // const G4double dydx[],
                          // const G4double Step,
                          G4double tau,
                          G4double yOut[] );
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 4; }

    void GetLastDydx( G4double dyDxLast[] );

    void   PrepareConstants();  // Initialise the values of the bi[][] array
   
  private :
    
    BogackiShampine45(const BogackiShampine45&);
    BogackiShampine45& operator=(const BogackiShampine45&);
    
    G4double *ak2, *ak3, *ak4,   *ak5, *ak6, *ak7, *ak8,
             *ak9, *ak10, *ak11, *yTemp,     *yIn;

    G4double *p[6];

    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector, *fLastDyDx,
             *fMidVector, *fMidError;
    // for DistChord calculations
    
    BogackiShampine45* fAuxStepper;  // For chord - until interpolation is proven
    bool   fPreparedInterpolation;

    // Class constants
    static bool     fPreparedConstants;   
    static G4double bi[12][7];
};

#endif /* defined(__Geant4__BogackiShampine45__) */
