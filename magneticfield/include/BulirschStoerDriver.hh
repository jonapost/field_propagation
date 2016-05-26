// $Id: G4MagIntegratorDriver.hh 69699 2013-05-13 08:50:30Z gcosmo $
//
//
// class BulirschStoerDriver
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and insures that
// the error is within acceptable bounds.

// History:
// - Created. D.Sorokins.
// --------------------------------------------------------------------

#ifndef BulirschStoerDriver_HH
#define BulirschStoerDriver_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4Mag_EqRhs.hh"
#include "G4EquationOfMotion.hh"

class BulirschStoerDriver{
public:

    BulirschStoerDriver(G4EquationOfMotion* _equation,
                        G4double _hmin,
                        G4int _integComp = 6,
                        G4int _verb = 1);

    ~BulirschStoerDriver();

    BulirschStoerDriver(const BulirschStoerDriver&) = delete;

    BulirschStoerDriver& operator=(const BulirschStoerDriver&) = delete;


     G4bool  AccurateAdvance(G4FieldTrack&  track,
                             G4double stepLen,
                             G4double eps,            // Requested y_err/hstep
                             G4double beginStep);  // Suggested 1st interval
       // Above drivers for integrator (Runge-Kutta) with stepsize control.
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps.
       // On output ystart is replaced by value at end of interval.
       // The concept is similar to the odeint routine from NRC p.721-722.

     G4bool  QuickAdvance(      G4FieldTrack& y_val,      // INOUT
                          const G4double     dydx[],
                                G4double     hstep,       // IN
                                G4double&    missDist,
                                G4double&    dyerr )  ;
        // QuickAdvance just tries one Step - it does not ensure accuracy.

     void GetDerivatives(const G4FieldTrack& track, G4double dydx[] );


private:

     void  OneGoodStep(G4double  y[],
                       const G4double  dydx[],
                       G4double& curveLength,
                       G4double htry,
                       G4double  eps,
                       G4double& hdid,
                       G4double& hnext ) ;
        // This takes one Step that is as large as possible while
        // satisfying the accuracy criterion of:
        // yerr < eps * |y_end-y_start|

     void ModifiedMidPoint(G4double y[], G4double dydx[], G4int order);

private:

     G4EquationOfMotion* equation;
     G4int verb;
     G4int integComp;
     G4double hmin;
};


#endif
