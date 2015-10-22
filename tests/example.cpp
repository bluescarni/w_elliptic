#include <complex>
#include <iostream>

// Include the w_elliptic.hpp header.
#include "w_elliptic.hpp"

int main()
{
    // Initialise a class for the computation of the Weierstrassian functions
    // in double precision, with real invariants g2 = 1 and g3 = 2.
    w_elliptic::we<double> w(1,2);

    // Computation of P with real argument.
    std::cout << w.P(1.2) << '\n';
    // Output: 0.923915...

    // Computation of P with complex argument.
    std::cout << w.P(std::complex<double>(1.2,3.4)) << '\n';
    // Output: (-1.07706...+0.268321...j)

    // Computation of zeta with real argument.
    std::cout << w.zeta(0.12) << '\n';
    // Output: 8.33330...

    // Computation of inverse P: will yield a pair of values in the
    // fundamental parallelogram.
    auto Pinv = w.Pinv(-4.);
    std::cout << Pinv[0] << '\n';
    // Output: (2.61673...+0.500504...j)
    std::cout << Pinv[1] << '\n';
    // Output: (1.30836...+1.95823...j)

    // Print a human-readable summary of the 'w' object:
    std::cout << w << '\n';
    // Output:
    // Invariants: [1,2]
    // Delta: -1712
    // Roots: [(0.898160...,0),(-0.449080...,0.595835...),(-0.449080...,-0.595835...)]
    // Periods: [(2.61673...,0),(1.30836...,2.45874...)]
    // etas: [(0.669458...,0),(0.334729...,-0.571542...)]
    // q: (0,0.0522398...)
}
