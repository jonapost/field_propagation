//
// class BulirschStoerDriver
//
// Class description:
//
// This is a driver class. It uses Bulirsch-Stoer algorithm to integrate ODE.
//
// History:
// - Created. D.Sorokin.
// --------------------------------------------------------------------

#ifndef BulirschStoerDriver_HH
#define BulirschStoerDriver_HH

#include "BaseDriver.hh"
#include "BSStepper.hh"

#include "ModifiedMidpoint.hh"
#include <boost/numeric/odeint.hpp>

class BulirschStoerDriver: public BaseDriver{
public:

    BulirschStoerDriver(G4double hminimum,
                        G4EquationOfMotion* equation,
                        G4int integratedComponents = 6,
                        G4int verb = 1);

    ~BulirschStoerDriver();

    BulirschStoerDriver(const BulirschStoerDriver&) = delete;

    BulirschStoerDriver& operator=(const BulirschStoerDriver&) = delete;


     virtual G4bool  AccurateAdvance(G4FieldTrack&  track,
                             G4double stepLen,
                             G4double eps,            // Requested y_err/hstep
                             G4double beginStep = 0);  // Suggested 1st interval
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps.
       // On output ystart is replaced by value at end of interval.


     virtual G4bool  QuickAdvance(      G4FieldTrack& y_val,      // INOUT
                          const G4double     dydx[],
                                G4double     hstep,
                                G4double&    missDist,
                                G4double&    dyerr )  ;
        // QuickAdvance just tries one Step - it does not ensure accuracy.

    virtual void GetDerivatives(const G4FieldTrack& track, G4double dydx[] );

    virtual G4MagIntegratorStepper* GetStepper();


private:

     void  OneGoodStep(G4double  y[],
                       const G4double  dydx[],
                       G4double& curveLength,
                       G4double htry,
                       G4double  eps,
                       G4double& hdid,
                       G4double& hnext) ;
        // This takes one Step that is as large as possible while
        // satisfying the accuracy criterion of:
        // yerr < eps * |y_end-y_start|

     //this is a dummy stepper to glue things up
     BSStepper* dummyStepper;
     ModifiedMidpoint mp;
     //boost::numeric::odeint::bulirsch_stoer<state_type> BulirschStoer;
};


#endif
