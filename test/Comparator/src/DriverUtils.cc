#include "DriverUtils.hh"

#include "templates.hh"

G4double relativeError(const G4double y[],
                       const G4double yerr[],
                       const G4double hstep)
{
    // square of displacement error
    G4double positionError2 =  sqr(yerr[0]) + sqr(yerr[1]) + sqr(yerr[2]);
    positionError2 /= sqr(hstep);

    // square of momentum vector difference
    G4double momentumError2 =  sqr(yerr[3]) + sqr(yerr[4]) + sqr(yerr[5]);
    G4double momentum2 =  sqr(y[3]) + sqr(y[4]) + sqr(y[5]);
    momentumError2 /= momentum2;

    // Square of maximum error
    return sqrt(std::max(positionError2, momentumError2));
}
