#ifndef UTILS_HH
#define UTILS_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"

enum class Type {
    Position = 0,
    Momentum = 3,
    Spin = 9
};

inline G4double extractValue2(const G4double array[], const Type& type)
{
    const size_t begin = static_cast<size_t>(type);
    return sqr(array[begin]) + sqr(array[begin+1]) + sqr(array[begin+2]);
}

inline G4ThreeVector makeVector(const G4double array[], const Type& type)
{
    const size_t begin = static_cast<size_t>(type);
    return G4ThreeVector(array[begin], array[begin + 1], array[begin + 2]);
}

#endif
