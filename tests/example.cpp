/***************************************************************************
 *   Copyright (C) 2015 by Francesco Biscani                               *
 *   bluescarni@gmail.com                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <complex>
#include <iostream>

#include "w_elliptic.hpp"

int main()
{
    // Initialise a class for the computation of the Weierstrassian functions
    // in double precision, with real invariants g2 = 1 and g3 = 2.
    w_elliptic::we<double> w(1,2);

    // Computation of P with real argument.
    std::cout << w.P(1.2) << '\n'; // Approximately 0.923915...

    // Computation of P with complex argument.
    std::cout << w.P(std::complex<double>(1.2,3.4)) << '\n'; // Approximately (-1.07706...+0.268321...j)

    // Computation of zeta with real argument.
    std::cout << w.zeta(0.12) << '\n'; // Approximately 8.33330...

    // Computation of inverse P: will yield a pair of values in the
    // fundamental parallelogram.
    auto Pinv = w.Pinv(-4.);
    std::cout << Pinv[0] << '\n'; // Approximately (2.61673...+0.500504...j)
    std::cout << Pinv[1] << '\n'; // Approximately (1.30836...+1.95823...j)
}
