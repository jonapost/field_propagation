// class G4VIntegrationDriver implementation by Dmitry Sorokin
// This is a base class for a driver algorithm.
//
// Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2016
//
//
// This code is made available subject to the Geant4 license, a copy of
//  which is available at http://www.geant4.org/geant4/license/
//
//  History
// -----------------------------
//  Created by Dmitry Sorokin 2016
//
//
///////////////////////////////////////////////////////////////////////////////

#ifndef G4VIntegrationDriver_HH
#define G4VIntegrationDriver_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4Mag_EqRhs.hh"
#include "G4EquationOfMotion.hh"
#include "G4MagIntegratorStepper.hh"

typedef std::pair<G4double, G4double> interpolationInterval;


class G4VIntegrationDriver{
public:

    //constructs driver without stepper
    G4VIntegrationDriver(G4double hminimum,
               G4EquationOfMotion* equation,
               G4int numberOfComponents = 6,
               G4int VerboseLevel = 1);

    //contructs driver with stepper
    G4VIntegrationDriver(G4double hminimum,
               G4MagIntegratorStepper *pItsStepper,
               G4int  numberOfComponents = 6,
               G4int  VerboseLevel = 1);

    virtual ~G4VIntegrationDriver();

    G4VIntegrationDriver(const G4VIntegrationDriver&) = delete;
    G4VIntegrationDriver& operator=(const G4VIntegrationDriver&) = delete;

    virtual G4bool  AccurateAdvance(G4FieldTrack&  track, G4double stepLen,
                                    G4double eps, G4double beginStep = 0) = 0;
    // Integrates ODE starting values y_current
    // from current s (s=s0) to s=s0+h with accuracy eps.
    // On output ystart is replaced by value at end of interval.


    virtual G4bool  QuickAdvance(G4FieldTrack& track, const G4double dydx[],
                                 G4double hstep, G4double& missDist,
                                 G4double& dyerr) = 0;
    // QuickAdvance just tries one Step - it does not ensure accuracy.

    virtual void OneGoodStep(G4double  ystart[],  const G4double  dydx[],
                             G4double& x, G4double htry,
                             G4double  eps, G4double& hdid,
                             G4double& hnext) = 0;
    // This takes one Step that is as large as possible while
    // satisfying the accuracy criterion of:
    // yerr < eps * |y_end-y_start|

    virtual G4double ComputeNewStepSize(G4double  errMaxNorm,    // normalised error
                                        G4double  hstepCurrent) = 0; // current step size
    // Taking the last step's normalised error, calculate
    // a step size for the next step.
    // Do not limit the next step's size within a factor of the
    // current one.

    virtual G4bool isDense() const = 0;

    virtual void DoStep(G4FieldTrack& track, G4double hstep, G4double eps);
    // Does one step with error control
    virtual void DoInterpolation(G4FieldTrack& track, G4double hstep, G4double eps = 0);
    virtual void Reset();


    inline void SetMinimumStep(G4double MinimumStep);
    inline void SetVerboseLevel(G4int VerboseLevel);
    inline void SetEquationOfMotion(G4EquationOfMotion* pEquation);
    inline void SetStepper(G4MagIntegratorStepper* pStepper);
    inline void SetInterpolationInterval(const interpolationInterval& interval);

    inline G4EquationOfMotion* GetEquationOfMotion() const;
    inline const G4MagIntegratorStepper* GetStepper() const;
    inline G4MagIntegratorStepper* GetStepper();
    inline G4double GetMinimumStep() const;
    inline G4int GetNumberOfVariables() const;
    inline G4double GetVerboseLevel() const;
    inline void GetDerivatives(const G4FieldTrack& track, G4double dydx[]) const;
    inline interpolationInterval& GetInterpolationInterval();

private:

     G4MagIntegratorStepper* fpStepper;
     G4EquationOfMotion* fpEquation;
     G4double fMinimumStep;
     G4int fNoIntegrationVariables;
     G4int fVerboseLevel;

     interpolationInterval fInterpolationInterval;
};

#include "G4VIntegrationDriver.icc"

#endif
