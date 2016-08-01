//
// class G4RKChordFinder
//
// Class description:
//
// Specialisation of G4VRevisedChordFinder class for
// Runge-Kutta drivers (G4MagInt_Driver)
//
// History:
// - Created: D. Sorokin
// --------------------------------------------------------------------

#ifndef G4RKChordFinder_HH
#define G4RKChordFinder_HH

#include "G4VRevisedChordFinder.hh"

class G4RKChordFinder : public G4VRevisedChordFinder
{
public:

    G4RKChordFinder(G4VIntegrationDriver* pIntegrationDriver ,
                    G4int staticsticsVerbosity = 1);

    //Constructor that creates defaults for all "children" classes.
    G4RKChordFinder(G4MagneticField* magField,
                    G4double stepMinimum = 1.0e-2, // * mm
                    G4MagIntegratorStepper* pItsStepper = nullptr,
                    G4int statisticsVerbosity = 1);

    ~G4RKChordFinder();

    G4RKChordFinder(const G4RKChordFinder& other);

    const G4RKChordFinder& operator = (const G4RKChordFinder& other);

    virtual G4double AdvanceChordLimited(G4FieldTrack& trackCurrent,
                                         G4double stepLen,
                                         G4double eps) override final;

    virtual void reset() override final;

    virtual G4bool DoStepForIntersection(G4FieldTrack&  track,
                                         G4double stepLen,
                                         G4double eps) override final;

    inline void SetEquationOfMotion(G4EquationOfMotion* newEquation);
    inline G4EquationOfMotion* GetEquationOfMotion();

private:

    G4double FindNextChord(const G4FieldTrack& yStart,
                           G4double stepMax,
                           G4FieldTrack& yEnd, // Endpoint
                           G4double& dyErrPos, // Error of endpoint
                           G4double epsStep,
                           G4double* pStepForAccuracy);

private:

    G4EquationOfMotion* fEquation;
    G4MagIntegratorStepper* fStepper;

    G4bool fAllocatedStepper;
    G4bool fAllocatedEquation;
};

#include "G4RKChordFinder.icc"

#endif
