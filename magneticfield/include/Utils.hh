#ifndef UTILS_HH
#define UTILS_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"

namespace magneticfield {

enum class Value3D {
    Position = 0,
    Momentum = 3,
    Spin = 9
};

enum class Value1D {
    KineticEnergy = 6,
    LabTime = 7,
    ProperTime = 8
};

inline G4double extractValue(const G4double array[], const Value1D& value);
inline G4double extractValue2(const G4double array[], const Value1D& value);

inline G4double extractValue(const G4double array[], const Value3D& value);
inline G4double extractValue2(const G4double array[], const Value3D& value);

inline G4ThreeVector makeVector(const G4double array[], const Value3D& value);

inline void copyValue(G4double dest[], const G4double src[], const Value1D& value);
inline void copyValue(G4double dest[], const G4double src[], const Value3D& value);

G4double relativeError(
    const G4double y[],
    const G4double yerr[],
    const G4double hstep,
    const G4double errorTolerance = 1);

} //magneticfield

#include "Utils.icc"

#endif
