#ifndef BETTIS_FUNCTIONS_HH
#define BETTIS_FUNCTIONS_HH

#include <cmath>

namespace magneticfield {
namespace internal {

template <int val>
struct factorial {
    enum {value = val * factorial<val-1>::value};
};

template <>
struct factorial<0> {
    enum {value = 1};
};

template <int order>
double reqursiveBettis(double nu)
{
    constexpr int j = order - 2;
    return (1. / factorial<j>::value - reqursiveBettis<j>(nu)) / (nu * nu);
}

template <>
double reqursiveBettis<0>(double nu)
{
    return cos(nu);
}

template <>
double reqursiveBettis<1>(double nu)
{
    return sin(nu) / nu;
}

template <int order>
double taylorBettis(double nu)
{
    static int MAX_ITER = 3;

    double result = 0;

    int sign = 1;
    double nuPow2k = 1;
    int jPlus2K = order;
    double jPlus2KFactorial = factorial<order>::value;

    for (int k = 0; k < MAX_ITER - 1; ++k) {
        result += sign * nuPow2k / jPlus2KFactorial;
        sign = -sign;
        nuPow2k *= nu * nu;
        jPlus2KFactorial *= (jPlus2K + 1) * (jPlus2K + 2);
        jPlus2K += 2;
    }
    result += sign * nuPow2k / jPlus2KFactorial;

    return result;
}

} // internal

template <int order>
double bettisFunction(double nu)
{
    return nu < 0.01 ?
        internal::taylorBettis<order>(nu) : internal::reqursiveBettis<order>(nu);
}

} // magneticfield

#endif
