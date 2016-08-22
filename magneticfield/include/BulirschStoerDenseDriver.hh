// BulirschStoerDenseDriver implementation by Dmitry Sorokin
// This driver class uses Bulirsch-Stoer (dense) method to integrate
// the equation of motion
//
// Driver with step size and order control as well as dense output.
// Very good if high precision and dense output is required.
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

#ifndef BulirschStoerDenseDriver_HH
#define BulirschStoerDenseDriver_HH

#include "G4VIntegrationDriver.hh"
#include "BulirschStoerDenseOut.hh"
#define USE_BOOST
#ifdef USE_BOOST
#include "boost/numeric/odeint.hpp"
#include <array>

typedef std::array<G4double, G4FieldTrack::ncompSVEC> state_type;
#endif

class BulirschStoerDenseDriver: public G4VIntegrationDriver{

public:
    BulirschStoerDenseDriver(G4double hminimum, G4EquationOfMotion* pequation,
                             G4int numberOfComponents = 6, G4int VerboseLevel = 1);

    ~BulirschStoerDenseDriver();

    BulirschStoerDenseDriver(const BulirschStoerDenseDriver&) = delete;
    BulirschStoerDenseDriver& operator=(const BulirschStoerDenseDriver&) = delete;


    virtual G4bool  AccurateAdvance(G4FieldTrack&  track,G4double stepLen,
                                    G4double eps,G4double hinitial = 0) override final;

    virtual G4bool QuickAdvance(G4FieldTrack& track, const G4double dydx[],
                                G4double hstep,G4double& missDist, 
                                G4double& dyerr ) override final;

    virtual void  OneGoodStep(G4double  y[], const G4double  dydx[],
                              G4double& curveLength, G4double htry,
                              G4double  eps, G4double& hdid,
                              G4double& hnext) override final;

    virtual G4double ComputeNewStepSize(G4double /*dyErr_relative*/, G4double lastStepLength ) override final;


    //dense output methods
    virtual G4bool isDense() const;
    virtual void DoStep(G4FieldTrack& track, G4double hstep, G4double eps) override final;
    virtual void DoInterpolation(G4FieldTrack& track, G4double hstep, G4double eps = 0) override final;

private:

    ModifiedMidpointDenseOut denseMidpoint;
    BulirschStoerDenseOut bulirschStoer;
#ifdef USE_BOOST
    boost::numeric::odeint::bulirsch_stoer_dense_out<state_type> boost_bulirsch_stoer;
#endif

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

    //for interpolation
    G4double eps_prev;
    G4double fNextStepSize;

};


#endif
