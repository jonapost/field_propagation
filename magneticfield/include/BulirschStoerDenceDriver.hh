// $Id: G4MagIntegratorDriver.hh 69699 2013-05-13 08:50:30Z gcosmo $
//
//
// class BulirschStoerDenceDriver
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and insures that
// the error is within acceptable bounds.

// History:
// - Created. D.Sorokins.
// --------------------------------------------------------------------

#ifndef BulirschStoerDenceDriver_HH
#define BulirschStoerDenceDriver_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4Mag_EqRhs.hh"
#include "G4EquationOfMotion.hh"

#include <functional>
#include <boost/array.hpp>

typedef boost::array<G4double,6> state_type;

class BulirschStoerDenceDriver{
public:

    BulirschStoerDenceDriver(G4EquationOfMotion* _equation,
                        G4int _integComp = 6,
                        G4int _verb = 1);

    ~BulirschStoerDenceDriver();

    BulirschStoerDenceDriver(const BulirschStoerDenceDriver&) = delete;

    BulirschStoerDenceDriver& operator=(const BulirschStoerDenceDriver&) = delete;


     G4bool  AccurateAdvance(G4FieldTrack&  track,
                             G4double stepLen,
                             G4double eps);  // Suggested 1st interval
       // Above drivers for integrator (Runge-Kutta) with stepsize control.
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps.
       // On output ystart is replaced by value at end of interval.
       // The concept is similar to the odeint routine from NRC p.721-722.

     G4double QuickAdvance(G4FieldTrack& y_val,
                                G4double     hstep,       // IN
                                G4double&    missDist,
                                G4double     epsStep )  ;
        // QuickAdvance just tries one Step - it does not ensure accuracy

     void GetDerivatives(const G4FieldTrack& track, G4double dydx[] );

private:

     const G4EquationOfMotion* equation;
     G4int verb;
     G4int integComp;
     G4double hmin;
     const G4double quickEps = 0.01;
     std::function<void(const state_type&, state_type&, double)> system;
};


#endif
