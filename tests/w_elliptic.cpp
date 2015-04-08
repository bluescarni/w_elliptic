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
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>
#include <fstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include <cmath>
#include <complex>

// TODO uniform the usage of these functions.
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

// Read a data file line by line, each line consisting of a set of comma-separated records.
static std::vector<std::vector<std::string>> read_data(const std::string &filename)
{
    std::vector<std::vector<std::string>> retval;
    std::ifstream infile(filename);
    std::string line;
    while (std::getline(infile, line)) {
        // Split the string.
        std::vector<std::string> tmp;
        boost::algorithm::split(tmp,line,boost::algorithm::is_any_of(","));
        retval.push_back(std::move(tmp));
    }
    // If the file was not found, throw.
    if (retval.empty()) {
        throw std::runtime_error("Loading of test data failed: make sure to run the test from the build/ directory.");
    }
    return retval;
}

static std::vector<std::vector<std::string>> test_01_data = read_data("../tests/test_01_data.txt");
static std::vector<std::vector<std::string>> test_02_data = read_data("../tests/test_02_data.txt");

typedef boost::mpl::vector<float,double,long double> real_types;

using size_type = std::vector<std::vector<std::string>>::size_type;

using namespace w_elliptic;

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
        complex_type eta;
        real_type max_eta_err(0), acc_eta_err(0);
        size_type max_eta_err_idx = 0u;
        real_type max_root_err(0), acc_root_err(0);
        size_type max_root_err_idx = 0u;
        real_type max_period_err(0), acc_period_err(0);
        size_type max_period_err_idx = 0u;
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
            complex_from_str(eta,v[12],v[13]);
            // Build the W object.
            we<real_type> w(g2,g3);
            // Compute the errors on the roots.
            for (std::size_t j = 0u; j < 3u; ++j) {
                if (roots[j] == real_type(0)) {
                    err = std::abs(w.roots()[j]);
                } else {
                    err = std::abs((roots[j]-w.roots()[j])/roots[j]);
                }
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
            // Errors on eta.
            err = std::abs((w.eta()-eta)/eta);
            acc_eta_err += err;
            if (err > max_eta_err) {
                max_eta_err = err;
                max_eta_err_idx = i;
            }
        }
        std::cout << "\tMax root error: " << max_root_err << " @ [g2=" << test_01_data[max_root_err_idx][0u] << ",g3=" << test_01_data[max_root_err_idx][1u] << "]\n";
        std::cout << "\tAverage root error: " << acc_root_err / real_type(test_01_data.size()) << '\n';
        std::cout << "\tMax period error: " << max_period_err << " @ [g2=" << test_01_data[max_period_err_idx][0u] << ",g3=" << test_01_data[max_period_err_idx][1u] << "]\n";
        std::cout << "\tAverage period error: " << acc_period_err / real_type(test_01_data.size()) << '\n';
        std::cout << "\tMax eta error: " << max_eta_err << " @ [g2=" << test_01_data[max_eta_err_idx][0u] << ",g3=" << test_01_data[max_eta_err_idx][1u] << "]\n";
        std::cout << "\tAverage eta error: " << acc_eta_err / real_type(test_01_data.size()) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_01)
{
    std::cout << "Testing the computation of roots and periods\n";
    std::cout << "============================================\n";
    boost::mpl::for_each<real_types>(tester_01());
}

struct tester_02
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        real_type g2, g3, x, P, P_comp, max_P_err = 0, acc_P_err = 0;
        size_type max_err_idx = 0;
        std::string max_err_x, max_err_P;
        for (size_type i = 0u; i < test_02_data.size(); ++i) {
            const auto &v = test_02_data[i];
            // Read the invariant values from mpmath.
            real_from_str(g2,v[0]);
            real_from_str(g3,v[1]);
            // Build the W object.
            we<real_type> w(g2,g3);
            for (decltype(v.size()) j = 2u; j < v.size(); j += 2u) {
                real_from_str(x,v[j]);
                P = w.P(x);
                real_from_str(P_comp,v[j + 1u]);
                if (std::abs((P-P_comp)/P_comp) > max_P_err) {
                    max_P_err = std::abs((P-P_comp)/P_comp);
                    max_err_idx = i;
                    max_err_x = v[j];
                    max_err_P = v[j+1];
                }
                acc_P_err += std::abs((P-P_comp)/P_comp);
            }
        }
        std::cout << "\tMax P error: " << max_P_err << " @ [g2=" << test_02_data[max_err_idx][0u] << ",g3=" << test_02_data[max_err_idx][1u] << ",x=" << max_err_x << ",P=" << max_err_P << "]\n";
        std::cout << "\tAverage P error: " << acc_P_err / (real_type(test_02_data.size())*100) << '\n';
//         std::cout << "\tMax period error: " << max_period_err << " @ [g2=" << test_01_data[max_period_err_idx][0u] << ",g3=" << test_01_data[max_period_err_idx][1u] << "]\n";
//         std::cout << "\tAverage period error: " << acc_period_err / real_type(test_01_data.size()) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_02)
{
    std::cout << "Testing the computation of real P\n";
    std::cout << "=================================\n";
    boost::mpl::for_each<real_types>(tester_02());
}
