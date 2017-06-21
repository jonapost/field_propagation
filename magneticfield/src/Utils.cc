#include "Utils.hh"

namespace {

template<class T>
size_t getFirstIndex(const T& value) {
    return static_cast<size_t>(value);
}

} //namespace

namespace magneticfield {

G4double extractValue(const G4double array[], const Value1D& value)
{
    const auto begin = getFirstIndex(value);
    return array[begin];
}

G4double extractValue2(const G4double array[], const Value1D& value)
{
    return sqr(extractValue(array, value));
}

G4double extractValue(const G4double array[], const Value3D &value)
{
    return sqrt(extractValue2(array, value));
}

G4double extractValue2(const G4double array[], const Value3D& value)
{
    const auto begin = getFirstIndex(value);
    return sqr(array[begin]) + sqr(array[begin+1]) + sqr(array[begin+2]);
}

G4ThreeVector makeVector(const G4double array[], const Value3D& value)
{
    const auto begin = getFirstIndex(value);
    return G4ThreeVector(array[begin], array[begin + 1], array[begin + 2]);
}

} //magneticfield
