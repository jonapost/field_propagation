//
// class BulirschStoerDenseDriver
//
// Class description:
//
// This is a driver class. It uses Bulirsch-Stoer algorithm to integrate ODE.
//
// History:
// - Created. D.Sorokin.
// --------------------------------------------------------------------

#ifndef BulirschStoerDenseDriver_HH
#define BulirschStoerDenseDriver_HH

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include "G4Mag_EqRhs.hh"
#include "G4EquationOfMotion.hh"

#include <functional>
#include <boost/array.hpp>

typedef boost::array<G4double,6> state_type;

class BulirschStoerDenseDriver{
public:
    BulirschStoerDenseDriver(G4EquationOfMotion* pequation,
                             G4int numberOfComponents = 6,
                             G4int statisticsVerbosity = 1);

    ~BulirschStoerDenseDriver();

    BulirschStoerDenseDriver(const BulirschStoerDenseDriver&) = delete;

    BulirschStoerDenseDriver& operator=(const BulirschStoerDenseDriver&) = delete;


    //integrates ODE
    void  AccurateAdvance(G4FieldTrack&  track,
                            G4double hstep,
                            G4double eps); //relative error

    //tries one Step with lower accuracy returns step did
    G4double QuickAdvance(G4FieldTrack& track,
                          G4double hstep,
                          G4double& chord, //chord between start and end points
                          G4double eps);

    G4double do_step(G4FieldTrack&  track,
                     G4double hstep,
                     G4double eps,
                     G4double dChordStep);

private:

    //dependent object
    G4EquationOfMotion* equation;

    //Number of Variables in integration
    const G4int nvar;

    //verbose level
    G4int verb;

    //epsilon for QuickAdvance
    const G4double quickEps;

    //system func for boost::odeint Bulirsch-Stoer solver
    std::function<void(const state_type&, state_type&, double)> system;
};


#endif
