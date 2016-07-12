//
// class BogackiShampine45
//
// Class description:
//
//  Implementation of the  BogackiShampine45 method
//  based on method 2 from RKSUITE
//
//
// History:
// - Created. D.Sorokin.
// --------------------------------------------------------------------


#ifndef BogackiShampine45_HH
#define BogackiShampine45_HH

#include "G4MagIntegratorStepper.hh"
#include "G4FieldTrack.hh"

class BogackiShampine45 : public G4MagIntegratorStepper
{
public:
    BogackiShampine45(G4EquationOfMotion *EqRhs,
                      G4int numberOfVariables = 6);

    ~BogackiShampine45();

    BogackiShampine45(const BogackiShampine45&) = delete;
    BogackiShampine45& operator=(const BogackiShampine45&) = delete;

    void Stepper(const G4double yInput[],
                  const G4double dydx[],
                        G4double h,
                        G4double yOutput[],
                        G4double yerr[] ) ;

    void Interpolate(G4double tau,G4double yOutput[]) const;

    G4double  DistChord()   const;

    G4int IntegratorOrder() const { return 4; }

  private :

    void Init();

    void FormInterpolation();

    G4double ak2[G4FieldTrack::ncompSVEC],
             ak3[G4FieldTrack::ncompSVEC],
             ak4[G4FieldTrack::ncompSVEC],
             ak5[G4FieldTrack::ncompSVEC],
             ak6[G4FieldTrack::ncompSVEC],
             ak7[G4FieldTrack::ncompSVEC],
             ak8[G4FieldTrack::ncompSVEC],
             ak9[G4FieldTrack::ncompSVEC],
             ak10[G4FieldTrack::ncompSVEC],
             ak11[G4FieldTrack::ncompSVEC],
             yTemp[G4FieldTrack::ncompSVEC],
             yIn[G4FieldTrack::ncompSVEC],
             yOut[G4FieldTrack::ncompSVEC],
             dydxIn[G4FieldTrack::ncompSVEC],
             dydxOut[G4FieldTrack::ncompSVEC],
             yError[G4FieldTrack::ncompSVEC];

    //coefs from RKSUITE
    G4double a[12][11];
    G4double b[9];
    G4double r[12][7];
    G4double e[8];
    //G4double c[12];

    //for interpolation
    G4double p[5][G4FieldTrack::ncompSVEC];

    mutable G4double yMid[G4FieldTrack::ncompSVEC];

    G4double fLastStepLength;

    G4bool fInterpolationPrepared;
};

#endif
