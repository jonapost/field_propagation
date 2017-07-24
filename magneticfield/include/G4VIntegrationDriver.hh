#ifndef G4VINTEGRATION_DRIVER_HH
#define G4VINTEGRATION_DRIVER_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4EquationOfMotion.hh"

class G4VIntegrationDriver {
public:
    G4VIntegrationDriver() = default;
    virtual ~G4VIntegrationDriver() = default;

    G4VIntegrationDriver(const G4VIntegrationDriver&) = delete;
    const G4VIntegrationDriver& operator = (const G4VIntegrationDriver&) = delete;

    virtual G4bool QuickAdvance(G4FieldTrack& track,   // INOUT
                                const G4double dydx[],
                                G4double hstep,
                                G4double& dchord_step,
                                G4double& dyerr) = 0;

    virtual G4bool AccurateAdvance(G4FieldTrack& track,
                                   G4double hstep,
                                   G4double eps,               // Requested y_err/hstep
                                   G4double hinitial = 0) = 0; // Suggested 1st interval

    virtual void GetDerivatives(const G4FieldTrack& track,
                                G4double dydx[]) const = 0;

    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) = 0;
    virtual G4EquationOfMotion* GetEquationOfMotion() = 0;

    // Taking the last step's normalised error, calculate
    // a step size for the next step.
    // Do not limit the next step's size within a factor of the
    // current one.
    virtual G4double ComputeNewStepSize(G4double errMaxNorm,    // normalised error
                                        G4double hstepCurrent) = 0; // current step size

    virtual void SetVerboseLevel(G4int level) = 0;
    virtual G4int GetVerboseLevel() const = 0;
};

#endif
