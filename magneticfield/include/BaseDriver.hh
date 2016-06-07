//
// class BaseDriver
//
// Class description:
//
// It is base class for BulirschStoerDriver and
// BulirschStoerDenseDriver
//
// History:
// - Created: D. Sorokin
// --------------------------------------------------------------------

#ifndef BaseDriver_HH
#define BaseDriver_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4Mag_EqRhs.hh"
#include "G4EquationOfMotion.hh"

#include <array>
#include <functional>

typedef std::array<G4double,G4FieldTrack::ncompSVEC> state_type;

class BaseDriver{
public:

    BaseDriver(G4EquationOfMotion* equation,
                        G4int integratedComponents = 6,
                        G4int verb = 1);

    virtual ~BaseDriver();

    BaseDriver(const BaseDriver&) = delete;

    BaseDriver& operator=(const BaseDriver&) = delete;


     virtual G4bool  AccurateAdvance(G4FieldTrack&  track,
                             G4double stepLen,
                             G4double eps,
                             G4double beginStep = 0) = 0;  // Suggested 1st interval
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps.
       // On output ystart is replaced by value at end of interval.


     virtual G4bool  QuickAdvance(      G4FieldTrack& y_val,      // INOUT
                          const G4double     dydx[],
                                G4double     hstep,       // IN
                                G4double&    missDist,
                                G4double&    dyerr )  = 0 ;
        // QuickAdvance just tries one Step - it does not ensure accuracy.

     virtual void GetDerivatives(const G4FieldTrack& track, G4double dydx[] ) = 0;

     G4double ComputeNewStepSize(double dyErr_relative,
                                 double lastStepLength );


//private:

     G4EquationOfMotion* fequation;
     G4int fverb;
     G4int fnvar;

     std::function<void(const state_type&, state_type&, G4double)> system;
     std::function<void(const state_type&, state_type&)> System;
};


#endif
