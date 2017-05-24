#include "DriverUtils.hh"

#include "globals.hh"
#include "templates.hh"

namespace {
    enum class Type {
        Position = 0,
        Momentum = 3,
        Spin = 9
    };

    G4double extractValue2(const G4double array[], const Type& type)
    {
        const size_t begin = static_cast<size_t>(type);
        return sqr(array[begin]) + sqr(array[begin+1]) + sqr(array[begin+2]);
    }
} //namespace


G4double relativeError(const G4double y[],
                       const G4double yerr[],
                       const G4double hstep)
{
    // square of displacement error
    G4double positionError2 =  extractValue2(yerr, Type::Position);
    positionError2 /= sqr(hstep);

    // square of momentum vector difference
    G4double momentum2 = extractValue2(y, Type::Momentum);
    G4double momentumError2 = 0;
    if (momentum2 > 0) {
        momentumError2 = extractValue2(yerr, Type::Momentum);
        momentumError2 /= momentum2;
    } else {
        G4Exception("driverUtils::relativeError()", "GeomField0003",
                    JustWarning, "momentum is zero");
    }

    // square of spin vector difference
    G4double spin2 = extractValue2(y, Type::Spin);
    G4double spinError2 = 0;
    if (spin2 > 0) {
        spinError2 = extractValue2(yerr, Type::Spin);
        spinError2 /= spin2;
    }

    // Square of maximum error
    return sqrt(std::max({positionError2, momentumError2, spinError2}));
}

