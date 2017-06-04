#ifndef RK547FEq1_HH
#define RK547FEq1_HH

#include "G4MagIntegratorStepper.hh"
#include "G4FieldTrack.hh"

class RK547FEq1 : public G4MagIntegratorStepper
{
public:
    RK547FEq1 (G4EquationOfMotion *EqRhs, G4int integrationVariables = 6);

    void Stepper (const G4double yInput[],
                  const G4double dydx[],
                        G4double hstep,
                        G4double yOutput[],
                        G4double yError[]);

    RK547FEq1 (const RK547FEq1&) = delete;
    RK547FEq1& operator = (const RK547FEq1&) = delete;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 4; }


private:
    void saveStep(const G4double yInput[],
                  const G4double dydx[],
                  const G4double hstep,
                  const G4double yOutput[]);

    void makeStep(const G4double yInput[],
                  const G4double dydx[],
                  const G4double hstep,
                        G4double yOutput[],
                        G4double yError[]) const;

    mutable G4double ak2_[G4FieldTrack::ncompSVEC],
                     ak3_[G4FieldTrack::ncompSVEC],
                     ak4_[G4FieldTrack::ncompSVEC],
                     ak5_[G4FieldTrack::ncompSVEC],
                     ak6_[G4FieldTrack::ncompSVEC],
                     ak7_[G4FieldTrack::ncompSVEC],

                     yMid_ [G4FieldTrack::ncompSVEC],
                     yErr_[G4FieldTrack::ncompSVEC],
                     yTemp_[G4FieldTrack::ncompSVEC];

    G4double yIn_[G4FieldTrack::ncompSVEC],
             yOut_[G4FieldTrack::ncompSVEC],
             dydx_[G4FieldTrack::ncompSVEC],
             dydxNext_[G4FieldTrack::ncompSVEC];
    G4double hstep_;
};

#endif
