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

#include "BaseDriver.hh"


#define private public
//#include "BulirschStoerDenseOut.hh"
#include "boost/numeric/odeint.hpp"
#define private private

#include "functional"

class BulirschStoerDenseDriver: public BaseDriver{
public:
    BulirschStoerDenseDriver(G4EquationOfMotion* pequation,
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


private:
    //epsilon for QuickAdvance
    const G4double quickEps;
    G4double tBegin,tEnd;
    boost::numeric::odeint::bulirsch_stoer_dense_out<state_type> theStepper;
    //BulirschStoerDenseOut theStepper;
};


#endif
