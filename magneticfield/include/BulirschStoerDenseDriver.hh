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


#define private public
//#include "BulirschStoerDenseOut.hh"
#include "boost/numeric/odeint.hpp"
#undef private

#include <array>
#include <functional>
typedef std::array<G4double,G4FieldTrack::ncompSVEC> state_type;

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
                            G4double beginStep = 0);  // Suggested 1st interval

    //tries one Step with lower accuracy returns step did
    virtual G4bool QuickAdvance(G4FieldTrack& track,
                                const G4double dydx[],
                                G4double hstep,
                                G4double& missDist, //chord between start and end points
                                G4double& dyerr );

    virtual void GetDerivatives(const G4FieldTrack& track, G4double dydx[] );

    virtual G4MagIntegratorStepper* GetStepper();


private:
    BSStepper* dummyStepper;
    //epsilon for QuickAdvance
    const G4double quickEps;
    G4double tBegin,tEnd;
    boost::numeric::odeint::bulirsch_stoer_dense_out<state_type> theStepper;
    //BulirschStoerDenseOut theStepper;


    std::function<void(const state_type& y, state_type& dydx, G4double t)> system;
};


#endif
