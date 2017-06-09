//
// class BulirschStoerDenseDriver
//
// Class description:
//
// This is a driver class. It uses Bulirsch-Stoer algorithm with
// dense output to integrate ODE.
//
// History:
// - Created. D.Sorokin.
// --------------------------------------------------------------------

#ifndef BulirschStoerDenseDriver_HH
#define BulirschStoerDenseDriver_HH

#include "G4VIntegrationDriver.hh"
#include "BSStepper.hh"


#include "BulirschStoerDenseOut.hh"
//#include "boost/numeric/odeint.hpp"


//#include <array>
//#include <functional>
//typedef std::array<G4double,G4FieldTrack::ncompSVEC> state_type;

class BulirschStoerDenseDriver: public G4VIntegrationDriver{
public:
    BulirschStoerDenseDriver(G4double hminimum,
                             G4EquationOfMotion* pequation,
                             G4int numberOfComponents = 6,
                             G4int statisticsVerbosity = 1);

    ~BulirschStoerDenseDriver();

    BulirschStoerDenseDriver(const BulirschStoerDenseDriver&) = delete;

    BulirschStoerDenseDriver& operator=(const BulirschStoerDenseDriver&) = delete;


    //integrates ODE
    virtual G4bool  AccurateAdvance(G4FieldTrack&  track,
                            G4double stepLen,
                            G4double eps,
                            G4double hinitial = 0);  // Suggested 1st interval

    //tries one Step with lower accuracy returns step did
    virtual G4bool QuickAdvance(G4FieldTrack& track,
                                const G4double dydx[],
                                G4double hstep,
                                G4double& missDist, //chord between start and end points
                                G4double& dyerr );

    virtual void GetDerivatives(const G4FieldTrack& track, G4double dydx[] );

    virtual G4MagIntegratorStepper* GetStepper();

    void  OneGoodStep(G4double  y[],
                      const G4double  dydx[],
                      G4double& curveLength,
                      G4double htry,
                      G4double  eps,
                      G4double& hdid,
                      G4double& hnext) ;

    virtual G4double ComputeNewStepSize(double /*dyErr_relative*/,
                                        double lastStepLength )
    {
        return lastStepLength;
    }


private:
    BSStepper* dummyStepper;
    //boost::numeric::odeint::bulirsch_stoer_dense_out<state_type> quickStepper;
    //boost::numeric::odeint::bulirsch_stoer_dense_out<state_type> accurateStepper;


    //std::function<void(const state_type& y, state_type& dydx, G4double t)> system;

    //this is a dummy stepper to glue things up
    ModifiedMidpointDenseOut denseMidpoint;
    BulirschStoerDenseOut bulirschStoer;
    //boost::numeric::odeint::bulirsch_stoer<state_type> BulirschStoer;



    G4double yIn[G4FieldTrack::ncompSVEC],
             yMid[G4FieldTrack::ncompSVEC],
             yMid2[G4FieldTrack::ncompSVEC],
             yOut[G4FieldTrack::ncompSVEC],
             yOut2[G4FieldTrack::ncompSVEC],
             yError[G4FieldTrack::ncompSVEC];



    G4double dydxCurrent[G4FieldTrack::ncompSVEC];
    G4double yCurrent[G4FieldTrack::ncompSVEC];

    G4double derivs[2][6][G4FieldTrack::ncompSVEC];

    // G4double diffs[4][2][G4FieldTrack::ncompSVEC];

    G4double tBegin,tEnd;

};


#endif
