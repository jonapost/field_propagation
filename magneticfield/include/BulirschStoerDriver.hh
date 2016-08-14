//
// class BulirschStoerDriver
//
// Class description:
//
// This is a driver class. It uses Bulirsch-Stoer algorithm to integrate ODE.
//
// History:
// - Created. D.Sorokin.
// --------------------------------------------------------------------

#ifndef BulirschStoerDriver_HH
#define BulirschStoerDriver_HH

#include "G4VIntegrationDriver.hh"
#include "ModifiedMidpoint.hh"
#include "ModifiedMidpointDenseOut.hh"
#include "BulirschStoer.hh"

class BulirschStoerDriver: public G4VIntegrationDriver{
public:

    BulirschStoerDriver(G4double hminimum,
                        G4EquationOfMotion* equation,
                        G4int integratedComponents = 6,
                        G4int verb = 1);

    ~BulirschStoerDriver();

    BulirschStoerDriver(const BulirschStoerDriver&) = delete;
    BulirschStoerDriver& operator=(const BulirschStoerDriver&) = delete;

    virtual G4bool  AccurateAdvance(G4FieldTrack&  track, G4double stepLen,
                                    G4double eps, G4double beginStep = 0) override final;



    virtual G4bool  QuickAdvance(G4FieldTrack& y_val, const G4double     dydx[],
                                 G4double hstep, G4double& missDist,
                                 G4double& dyerr ) override final;

    virtual void  OneGoodStep(G4double  y[],  const G4double  dydx[],
                              G4double& curveLength, G4double htry,
                              G4double  eps, G4double& hdid,
                              G4double& hnext) override final;

    virtual G4double ComputeNewStepSize(G4double  errMaxNorm,    // normalised error
                                        G4double  hstepCurrent) override final; // current step size

    virtual G4bool isDense() const override final;

private:
     ModifiedMidpoint modifiedMidpoint;
     ModifiedMidpointDenseOut denseMidpoint;
     BulirschStoer bulirschStoer;


     G4double yIn[G4FieldTrack::ncompSVEC],
              yMid[G4FieldTrack::ncompSVEC],
              yMid2[G4FieldTrack::ncompSVEC],
              yOut[G4FieldTrack::ncompSVEC],
              yOut2[G4FieldTrack::ncompSVEC],
              yError[G4FieldTrack::ncompSVEC];


     G4double dydxCurrent[G4FieldTrack::ncompSVEC];
     G4double yCurrent[G4FieldTrack::ncompSVEC];

     G4double derivs[2][6][G4FieldTrack::ncompSVEC];

     G4double diffs[4][2][G4FieldTrack::ncompSVEC];


     const G4int interval_sequence[2];
     const G4double fcoeff;

};


#endif
