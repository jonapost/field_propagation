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

G4double relativeError(
    const G4double y[],
    const G4double yerr[],
    const G4double hstep,
    const G4double errorTolerance,
    G4IntegrationObserver& observer)
{
    const G4double errorTolerance2 = sqr(errorTolerance);
    const G4double hstep2 = sqr(hstep);

    // square of displacement error
    G4double positionError2 =  extractValue2(yerr, Value3D::Position);
    positionError2 /= hstep2 * errorTolerance2;

    // square of momentum vector difference
    const G4double momentum2 = extractValue2(y, Value3D::Momentum);
    G4double momentumError2 = 0;
    if (momentum2 > 0) {
        momentumError2 = extractValue2(yerr, Value3D::Momentum);
        momentumError2 /= momentum2 * errorTolerance2;
    } else {
        G4Exception("G4MagInt_Driver::relativeError()", "GeomField0003",
                    JustWarning, "momentum is zero");
    }

    // square of spin vector difference
    const G4double spin2 = extractValue2(y, Value3D::Spin);
    G4double spinError2 = 0;
    if (spin2 > 0) {
        spinError2 = extractValue2(yerr, Value3D::Spin);
        spinError2 /= spin2 * errorTolerance2;
    }

    observer.onRelativeError(hstep, positionError2, momentumError2 * hstep2);

    // Square of maximum error
    return sqrt(std::max({positionError2, momentumError2, spinError2}));
}

} //magneticfield
