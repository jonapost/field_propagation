// class G4BogackiShampine45DenseDriver implementation by Dmitry Sorokin
// Implementaton of dense driver for Bogacki-Shampine45 algorithm
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

#ifndef G4BogackiShampine45DenseDriver_HH
#define G4BogackiShampine45DenseDriver_HH

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

#endif
