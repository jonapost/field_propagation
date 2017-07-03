#ifndef RK547FEq2_HH
#define RK547FEq2_HH

#include "G4MagIntegratorStepper.hh"
#include "G4FieldTrack.hh"

class RK547FEq2 : public G4MagIntegratorStepper {
public:
    RK547FEq2(G4EquationOfMotion* EqRhs, G4int integrationVariables = 6);

    virtual void Stepper(
        const G4double yInput[],
        const G4double dydx[],
        G4double hstep,
        G4double yOutput[],
        G4double yError[]) override;

    void Stepper(const G4double yInput[],
        const G4double dydx[],
        G4double hstep,
        G4double yOutput[],
        G4double yError[],
        G4double dydxOutput[]);

    RK547FEq2 (const RK547FEq2&) = delete;
    RK547FEq2& operator = (const RK547FEq2&) = delete;

    virtual G4double DistChord() const override;
    virtual G4int IntegratorOrder() const override { return 4; }

private:
    void makeStep(
        const G4double yInput[],
        const G4double dydx[],
        const G4double hstep,
        G4double yOutput[],
        G4double* dydxOutput = nullptr,
        G4double* yError = nullptr) const;

    G4double fyIn[G4FieldTrack::ncompSVEC],
             fdydx[G4FieldTrack::ncompSVEC],
             fyOut[G4FieldTrack::ncompSVEC],
             fdydxOut[G4FieldTrack::ncompSVEC];
    G4double fhstep;
};

#endif
