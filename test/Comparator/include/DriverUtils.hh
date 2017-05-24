#ifndef DRIVERUTILS_HH
#define DRIVERUTILS_HH

#include "G4Types.hh"


G4double relativeError(const G4double y[],
                       const G4double yerr[],
                       const G4double hstep);

#endif
