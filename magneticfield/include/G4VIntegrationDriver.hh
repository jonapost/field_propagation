//
// class G4VIntegrationDriver
//
// Class description:
//
// It is a base driver class
//
// History:
// - Created: D. Sorokin
// --------------------------------------------------------------------

#ifndef G4VIntegrationDriver_HH
#define G4VIntegrationDriver_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4Mag_EqRhs.hh"
#include "G4EquationOfMotion.hh"
#include "G4MagIntegratorStepper.hh"


class G4VIntegrationDriver{
public:

    //constructor for BulirschStoer driver
    G4VIntegrationDriver(G4double hminimum,
               G4EquationOfMotion* equation,
               G4int numberOfComponents = 6,
               G4int statisticsVerbosity = 1);

    G4VIntegrationDriver(G4double hminimum,
               G4MagIntegratorStepper *pItsStepper,
               G4int  numberOfComponents = 6,
               G4int  statisticsVerbosity = 1);

    virtual ~G4VIntegrationDriver();

    G4VIntegrationDriver(const G4VIntegrationDriver&) = delete;

    G4VIntegrationDriver& operator=(const G4VIntegrationDriver&) = delete;


    virtual G4bool  AccurateAdvance(G4FieldTrack&  track,
                             G4double stepLen,
                             G4double eps,
                             G4double beginStep = 0) = 0;  // Suggested 1st interval
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps.
       // On output ystart is replaced by value at end of interval.


    virtual G4bool  QuickAdvance(      G4FieldTrack& y_val,
                          const G4double     dydx[],
                                G4double     hstep,
                                G4double&    missDist,
                                G4double&    dyerr )  = 0 ;
     // QuickAdvance just tries one Step - it does not ensure accuracy.

    virtual void GetDerivatives(const G4FieldTrack& track, G4double dydx[] ) = 0;

    virtual G4double ComputeNewStepSize(double dyErr_relative,
                                        double lastStepLength ) = 0;


    virtual G4MagIntegratorStepper* GetStepper() = 0;


    inline void SetMinStep(G4double newStep);
    inline void SetVerboseLevel(G4int newVerb);

    inline G4EquationOfMotion* GetEquationOfMotion();
    inline void SetEquationOfMotion(G4EquationOfMotion* newEquation);
    inline G4double GetMinStep() const;
    inline G4int GetNumberOfVariables() const;
    inline G4double GetVerboseLevel() const;

private:

     G4MagIntegratorStepper* fStepper;
     G4EquationOfMotion* fequation;
     G4double fMinimumStep;
     G4int fnvar;
     G4int fverb;
};

#include "G4VIntegrationDriver.icc"

#endif
