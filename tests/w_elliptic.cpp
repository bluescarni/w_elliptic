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

#include "../src/w_elliptic.hpp"

#define BOOST_TEST_MODULE w_elliptic_test
#include <boost/test/unit_test.hpp>

#include <array>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>
#include <string>
#include <typeinfo>
#include <vector>

#include <cmath>
#include <complex>

static inline void real_from_str(float &x, const std::string &s)
{
    try {
        x = std::stof(s);
    } catch (...) {
        x = 0.f;
    }
}

static inline void complex_from_str(std::complex<float> &c, const std::string &re, const std::string &im)
{
    float R, I;
    real_from_str(R,re);
    real_from_str(I,im);
    c = std::complex<float>(R,I);
}

static inline void real_from_str(double &x, const std::string &s)
{
    x = std::stod(s);
}

static inline void complex_from_str(std::complex<double> &c, const std::string &re, const std::string &im)
{
    double R = std::stod(re);
    double I = std::stod(im);
    c = std::complex<double>(R,I);
}

static inline void real_from_str(long double &x, const std::string &s)
{
    x = std::stold(s);
}

static inline void complex_from_str(std::complex<long double> &c, const std::string &re, const std::string &im)
{
    long double R = std::stold(re);
    long double I = std::stold(im);
    c = std::complex<long double>(R,I);
}

using namespace w_elliptic;

typedef boost::mpl::vector<float,double,long double> real_types;

using size_type = std::vector<std::vector<std::string>>::size_type;

extern std::vector<std::vector<std::string>> test_01_data;

struct tester_01
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        using complex_type = typename we<real_type>::complex_type;
        real_type g2, g3, err;
        std::array<complex_type,3> roots;
        std::array<complex_type,2> periods;
        real_type max_root_err(0), acc_root_err(0);
        std::size_t max_root_err_idx = 0u;
        real_type max_period_err(0), acc_period_err(0);
        std::size_t max_period_err_idx = 0u;
        for (size_type i = 0u; i < test_01_data.size(); ++i) {
            const auto &v = test_01_data[i];
            // Read the values from mpmath.
            real_from_str(g2,v[0]);
            real_from_str(g3,v[1]);
            complex_from_str(roots[0],v[2],v[3]);
            complex_from_str(roots[1],v[4],v[5]);
            complex_from_str(roots[2],v[6],v[7]);
            complex_from_str(periods[0],v[8],v[9]);
            complex_from_str(periods[1],v[10],v[11]);
            // Build the W object.
            we<real_type> w(g2,g3);
            // Compute the errors on the roots.
            for (std::size_t j = 0u; j < 3u; ++j) {
                err = std::abs((roots[j]-w.roots()[j])/roots[j]);
                acc_root_err += err;
                if (err > max_root_err) {
                    max_root_err_idx = i;
                    max_root_err = err;
                }
            }
            // Compute the errors on the periods.
            for (std::size_t j = 0u; j < 2u; ++j) {
                err = std::abs((periods[j]-w.periods()[j])/periods[j]);
                acc_period_err += err;
                if (err > max_period_err) {
                    max_period_err_idx = i;
                    max_period_err = err;
                }
            }
        }
        std::cout << "\tMax root error: " << max_root_err << " @ [g2=" << test_01_data[max_root_err_idx][0u] << ",g3=" << test_01_data[max_root_err_idx][1u] << "]\n";
        std::cout << "\tAverage root error: " << acc_root_err / real_type(test_01_data.size()) << '\n';
        std::cout << "\tMax period error: " << max_period_err << " @ [g2=" << test_01_data[max_period_err_idx][0u] << ",g3=" << test_01_data[max_period_err_idx][1u] << "]\n";
        std::cout << "\tAverage period error: " << acc_period_err / real_type(test_01_data.size()) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_01)
{
    std::cout << "Testing the computation of roots and periods\n";
    std::cout << "============================================\n";
    boost::mpl::for_each<real_types>(tester_01());
}
