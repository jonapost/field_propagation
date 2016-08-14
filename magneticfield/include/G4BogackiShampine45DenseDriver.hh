#include "G4MagIntegratorDriver.hh"

class G4BogackiShampine45DenseDriver : public G4MagInt_Driver
{
public:
    G4BogackiShampine45DenseDriver(G4double hminimum, G4EquationOfMotion* pequation,
                                   G4int numberOfComponents = 6, G4int VerboseLevel = 1);
    virtual ~G4BogackiShampine45DenseDriver();

    G4BogackiShampine45DenseDriver(const G4BogackiShampine45DenseDriver& ) = delete;
    const G4BogackiShampine45DenseDriver& operator = (const G4BogackiShampine45DenseDriver& ) = delete;

    virtual G4bool isDense() const override final;
    virtual void DoStep(G4FieldTrack& track, G4double hstep, G4double eps) override final;
    virtual void DoInterpolation(G4FieldTrack& track, G4double hstep, G4double eps) override final;

private:

    //for interpolation
    G4double eps_prev;
    G4double fNextStepSize;

    G4double yCurrent[G4FieldTrack::ncompSVEC],
             dydxCurrent[G4FieldTrack::ncompSVEC];

};
