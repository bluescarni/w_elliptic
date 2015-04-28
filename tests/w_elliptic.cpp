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

#include <algorithm>
#include <array>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>
#include <cmath>
#include <fstream>
#include <random>
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

// Test data - these are tables of real numbers in string form.
static std::vector<std::vector<std::string>> test_01_data = read_data("../tests/test_01_data.txt");
static std::vector<std::vector<std::string>> test_02_data = read_data("../tests/test_02_data.txt");
static std::vector<std::vector<std::string>> test_03_data = read_data("../tests/test_03_data.txt");
static std::vector<std::vector<std::string>> test_04_data = read_data("../tests/test_04_data.txt");
static std::vector<std::vector<std::string>> test_05_data = read_data("../tests/test_05_data.txt");
static std::vector<std::vector<std::string>> test_06_data = read_data("../tests/test_06_data.txt");
static std::vector<std::vector<std::string>> test_07_data = read_data("../tests/test_07_data.txt");
static std::vector<std::vector<std::string>> test_09_data = read_data("../tests/test_09_data.txt");
static std::vector<std::vector<std::string>> test_11_data = read_data("../tests/test_11_data.txt");

// RNG.
static std::mt19937 rng;

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
            // Order the roots according to the expected convention.
            auto w_roots = w.roots();
            if (w.Delta() < real_type(0)) {
                std::sort(w_roots.begin(),w_roots.end(),[](const complex_type &c1, const complex_type &c2) {
                    return c1.imag() > c2.imag();
                });
                // The second root is real.
                w_roots[1] = complex_type(w_roots[1].real(),real_type(0));
            } else {
                std::sort(w_roots.begin(),w_roots.end(),[](const complex_type &c1, const complex_type &c2) {
                    return c1.real() > c2.real();
                });
                // All the roots are real.
                std::for_each(w_roots.begin(),w_roots.end(),[](complex_type &c) {
                    c = complex_type(c.real(),real_type(0));
                });
            }
            // Compute the errors on the roots.
            for (std::size_t j = 0u; j < 3u; ++j) {
                if (roots[j] == real_type(0)) {
                    err = std::abs(w_roots[j]);
                } else {
                    err = std::abs((roots[j]-w_roots[j])/roots[j]);
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
            err = std::abs((w.etas()[0]-eta)/eta);
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
    std::cout << "Testing the computation of roots, periods and eta\n";
    std::cout << "=================================================\n";
    boost::mpl::for_each<real_types>(tester_01());
    std::cout << "\n\n\n";
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
    }
};

BOOST_AUTO_TEST_CASE(test_02)
{
    std::cout << "Testing the computation of real P\n";
    std::cout << "=================================\n";
    boost::mpl::for_each<real_types>(tester_02());
    std::cout << "\n\n\n";
}

struct tester_03
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        using complex_type = typename we<real_type>::complex_type;
        real_type g2, g3, max_P_err = 0, acc_P_err = 0;
        complex_type c, P, P_comp;
        size_type max_err_idx = 0;
        std::string max_err_x, max_err_P;
        for (size_type i = 0u; i < test_03_data.size(); ++i) {
            const auto &v = test_03_data[i];
            // Read the invariant values from mpmath.
            real_from_str(g2,v[0]);
            real_from_str(g3,v[1]);
            // Build the W object.
            we<real_type> w(g2,g3);
            for (decltype(v.size()) j = 2u; j < v.size(); j += 4u) {
                complex_from_str(c,v[j],v[j + 1u]);
                P = w.P(c);
                complex_from_str(P_comp,v[j + 2u],v[j + 3u]);
                if (std::abs((P-P_comp)/P_comp) > max_P_err) {
                    max_P_err = std::abs((P-P_comp)/P_comp);
                    max_err_idx = i;
                    max_err_x = std::string("(") + v[j] + "," + v[j + 1u] + ")";
                    max_err_P = std::string("(") + v[j + 2u] + "," + v[j + 3u] + ")";
                }
                acc_P_err += std::abs((P-P_comp)/P_comp);
            }
        }
        std::cout << "\tMax P error: " << max_P_err << " @ [g2=" << test_03_data[max_err_idx][0u] << ",g3=" << test_03_data[max_err_idx][1u] << ",c=" << max_err_x << ",P=" << max_err_P << "]\n";
        std::cout << "\tAverage P error: " << acc_P_err / (real_type(test_03_data.size())*100) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_03)
{
    std::cout << "Testing the computation of complex P\n";
    std::cout << "====================================\n";
    boost::mpl::for_each<real_types>(tester_03());
    std::cout << "\n\n\n";
}

struct tester_04
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        real_type g2, g3, x, P, P_comp, max_P_err = 0, acc_P_err = 0;
        size_type max_err_idx = 0;
        std::string max_err_x, max_err_P;
        for (size_type i = 0u; i < test_04_data.size(); ++i) {
            const auto &v = test_04_data[i];
            // Read the invariant values from mpmath.
            real_from_str(g2,v[0]);
            real_from_str(g3,v[1]);
            // Build the W object.
            we<real_type> w(g2,g3);
            for (decltype(v.size()) j = 2u; j < v.size(); j += 2u) {
                real_from_str(x,v[j]);
                P = w.Pprime(x);
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
        std::cout << "\tMax P' error: " << max_P_err << " @ [g2=" << test_04_data[max_err_idx][0u] << ",g3=" << test_04_data[max_err_idx][1u] << ",x=" << max_err_x << ",P'=" << max_err_P << "]\n";
        std::cout << "\tAverage P' error: " << acc_P_err / (real_type(test_04_data.size())*100) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_04)
{
    std::cout << "Testing the computation of real P'\n";
    std::cout << "==================================\n";
    boost::mpl::for_each<real_types>(tester_04());
    std::cout << "\n\n\n";
}

struct tester_05
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        using complex_type = typename we<real_type>::complex_type;
        real_type g2, g3, max_P_err = 0, acc_P_err = 0;
        complex_type c, P, P_comp;
        size_type max_err_idx = 0;
        std::string max_err_x, max_err_P;
        for (size_type i = 0u; i < test_05_data.size(); ++i) {
            const auto &v = test_05_data[i];
            // Read the invariant values from mpmath.
            real_from_str(g2,v[0]);
            real_from_str(g3,v[1]);
            // Build the W object.
            we<real_type> w(g2,g3);
            for (decltype(v.size()) j = 2u; j < v.size(); j += 4u) {
                complex_from_str(c,v[j],v[j + 1u]);
                P = w.Pprime(c);
                complex_from_str(P_comp,v[j + 2u],v[j + 3u]);
                if (std::abs((P-P_comp)/P_comp) > max_P_err) {
                    max_P_err = std::abs((P-P_comp)/P_comp);
                    max_err_idx = i;
                    max_err_x = std::string("(") + v[j] + "," + v[j + 1u] + ")";
                    max_err_P = std::string("(") + v[j + 2u] + "," + v[j + 3u] + ")";
                }
                acc_P_err += std::abs((P-P_comp)/P_comp);
            }
        }
        std::cout << "\tMax P' error: " << max_P_err << " @ [g2=" << test_05_data[max_err_idx][0u] << ",g3=" << test_05_data[max_err_idx][1u] << ",c=" << max_err_x << ",P'=" << max_err_P << "]\n";
        std::cout << "\tAverage P' error: " << acc_P_err / (real_type(test_05_data.size())*100) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_05)
{
    std::cout << "Testing the computation of complex P'\n";
    std::cout << "=====================================\n";
    boost::mpl::for_each<real_types>(tester_05());
    std::cout << "\n\n\n";
}

struct tester_06
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        real_type g2, g3, x, P, P_comp, max_P_err = 0, acc_P_err = 0;
        size_type max_err_idx = 0;
        std::string max_err_x, max_err_P;
        for (size_type i = 0u; i < test_06_data.size(); ++i) {
            const auto &v = test_06_data[i];
            // Read the invariant values from mpmath.
            real_from_str(g2,v[0]);
            real_from_str(g3,v[1]);
            // Build the W object.
            we<real_type> w(g2,g3);
            for (decltype(v.size()) j = 2u; j < v.size(); j += 2u) {
                real_from_str(x,v[j]);
                P = w.zeta(x);
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
        std::cout << "\tMax zeta error: " << max_P_err << " @ [g2=" << test_06_data[max_err_idx][0u] << ",g3=" << test_06_data[max_err_idx][1u] << ",x=" << max_err_x << ",zeta=" << max_err_P << "]\n";
        std::cout << "\tAverage zeta error: " << acc_P_err / (real_type(test_06_data.size())*100) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_06)
{
    std::cout << "Testing the computation of real zeta\n";
    std::cout << "====================================\n";
    boost::mpl::for_each<real_types>(tester_06());
    std::cout << "\n\n\n";
}

struct tester_07
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        using complex_type = typename we<real_type>::complex_type;
        real_type g2, g3, max_P_err = 0, acc_P_err = 0;
        complex_type c, P, P_comp;
        size_type max_err_idx = 0;
        std::string max_err_x, max_err_P;
        for (size_type i = 0u; i < test_07_data.size(); ++i) {
            const auto &v = test_07_data[i];
            // Read the invariant values from mpmath.
            real_from_str(g2,v[0]);
            real_from_str(g3,v[1]);
            // Build the W object.
            we<real_type> w(g2,g3);
            for (decltype(v.size()) j = 2u; j < v.size(); j += 4u) {
                complex_from_str(c,v[j],v[j + 1u]);
                P = w.zeta(c);
                complex_from_str(P_comp,v[j + 2u],v[j + 3u]);
                if (std::abs((P-P_comp)/P_comp) > max_P_err) {
                    max_P_err = std::abs((P-P_comp)/P_comp);
                    max_err_idx = i;
                    max_err_x = std::string("(") + v[j] + "," + v[j + 1u] + ")";
                    max_err_P = std::string("(") + v[j + 2u] + "," + v[j + 3u] + ")";
                }
                acc_P_err += std::abs((P-P_comp)/P_comp);
            }
        }
        std::cout << "\tMax zeta error: " << max_P_err << " @ [g2=" << test_07_data[max_err_idx][0u] << ",g3=" << test_07_data[max_err_idx][1u] << ",c=" << max_err_x << ",zeta=" << max_err_P << "]\n";
        std::cout << "\tAverage zeta error: " << acc_P_err / (real_type(test_07_data.size())*100) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_07)
{
    std::cout << "Testing the computation of complex zeta\n";
    std::cout << "=======================================\n";
    boost::mpl::for_each<real_types>(tester_07());
    std::cout << "\n\n\n";
}

struct tester_08
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        using complex_type = typename we<real_type>::complex_type;
        std::uniform_real_distribution<double> rdist(-100.,100.);
        std::uniform_real_distribution<double> g2dist(-10.,20.);
        std::uniform_real_distribution<double> g3dist(-10.,10.);
        real_type max_Pinv_err = 0., max_g2 = 0., max_g3 = 0., acc_err = 0.;
        complex_type max_P(real_type(0.),real_type(0.)), max_P_inv(max_P);
        unsigned long counter = 0;
        for (int i = 0; i < 20; ++i) {
            real_type g2(static_cast<real_type>(g2dist(rng))), g3(static_cast<real_type>(g3dist(rng)));
            we<real_type> w(g2,g3);
            for (int j = 0; j < 1000; ++j) {
                complex_type P(real_type(rdist(rng)),real_type(rdist(rng)));
                auto Pinv = w.Pinv(P);
                auto err = std::abs((P-w.P(Pinv))/P);
                if (err > max_Pinv_err) {
                    max_Pinv_err = err;
                    max_g2 = g2;
                    max_g3 = g3;
                    max_P = P;
                    max_P_inv = Pinv;
                }
                acc_err += err;
                ++counter;
            }
        }
        // Add a couple of tests with g2/g3 == 0.
        for (int i = 0; i < 4; ++i) {
            real_type g2(0), g3(static_cast<real_type>(g3dist(rng)));
            we<real_type> w(g2,g3);
            for (int j = 0; j < 1000; ++j) {
                complex_type P(real_type(rdist(rng)),real_type(rdist(rng)));
                auto Pinv = w.Pinv(P);
                auto err = std::abs((P-w.P(Pinv))/P);
                if (err > max_Pinv_err) {
                    max_Pinv_err = err;
                    max_g2 = g2;
                    max_g3 = g3;
                    max_P = P;
                    max_P_inv = Pinv;
                }
                acc_err += err;
                ++counter;
            }
        }
        for (int i = 0; i < 4; ++i) {
            real_type g2(static_cast<real_type>(g2dist(rng))), g3(0);
            we<real_type> w(g2,g3);
            for (int j = 0; j < 1000; ++j) {
                complex_type P(real_type(rdist(rng)),real_type(rdist(rng)));
                auto Pinv = w.Pinv(P);
                auto err = std::abs((P-w.P(Pinv))/P);
                if (err > max_Pinv_err) {
                    max_Pinv_err = err;
                    max_g2 = g2;
                    max_g3 = g3;
                    max_P = P;
                    max_P_inv = Pinv;
                }
                acc_err += err;
                ++counter;
            }
        }
        std::cout << "\tMax Pinv error: " << max_Pinv_err << " @ [g2=" << max_g2 << ",g3=" << max_g3 << ",P=" << max_P  << ",Pinv=" << max_P_inv << "]\n";
        std::cout << "\tAverage Pinv error: " << acc_err / real_type(counter) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_08)
{
    std::cout << "Testing the computation of inverse P\n";
    std::cout << "====================================\n";
    boost::mpl::for_each<real_types>(tester_08());
    std::cout << "\n\n\n";
}

struct tester_09
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        using complex_type = typename we<real_type>::complex_type;
        real_type g2, g3, max_P_err = 0, acc_P_err = 0;
        complex_type c, P, P_comp;
        size_type max_err_idx = 0;
        std::string max_err_x, max_err_P;
        for (size_type i = 0u; i < test_09_data.size(); ++i) {
            const auto &v = test_09_data[i];
            // Read the invariant values from mpmath.
            real_from_str(g2,v[0]);
            real_from_str(g3,v[1]);
            // Build the W object.
            we<real_type> w(g2,g3);
            for (decltype(v.size()) j = 2u; j < v.size(); j += 4u) {
                complex_from_str(c,v[j],v[j + 1u]);
                P = w.sigma(c);
                if (!detail::isfinite(P)) {
                    continue;
                }
                complex_from_str(P_comp,v[j + 2u],v[j + 3u]);
                if (std::abs((P-P_comp)/P_comp) > max_P_err) {
                    max_P_err = std::abs((P-P_comp)/P_comp);
                    max_err_idx = i;
                    max_err_x = std::string("(") + v[j] + "," + v[j + 1u] + ")";
                    max_err_P = std::string("(") + v[j + 2u] + "," + v[j + 3u] + ")";
                }
                acc_P_err += std::abs((P-P_comp)/P_comp);
            }
        }
        std::cout << "\tMax sigma error: " << max_P_err << " @ [g2=" << test_09_data[max_err_idx][0u] << ",g3=" << test_09_data[max_err_idx][1u] << ",c=" << max_err_x << ",sigma=" << max_err_P << "]\n";
        std::cout << "\tAverage sigma error: " << acc_P_err / (real_type(test_09_data.size())*100) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_09)
{
    std::cout << "Testing the computation of complex sigma\n";
    std::cout << "========================================\n";
    boost::mpl::for_each<real_types>(tester_09());
    std::cout << "\n\n\n";
}

struct tester_10
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        using complex_type = typename we<real_type>::complex_type;
        std::uniform_real_distribution<double> rdist(-2.,2.);
        std::uniform_real_distribution<double> g2dist(-10.,20.);
        std::uniform_real_distribution<double> g3dist(-10.,10.);
        real_type max_err = 0., max_g2 = 0., max_g3 = 0., acc_err = 0.;
        complex_type max_z(0);
        unsigned long counter = 0;
        for (int i = 0; i < 20; ++i) {
            real_type g2(static_cast<real_type>(g2dist(rng))), g3(static_cast<real_type>(g3dist(rng)));
            we<real_type> w(g2,g3);
            for (int j = 0; j < 1000; ++j) {
                complex_type z(real_type(rdist(rng))*w.periods()[0] + real_type(rdist(rng))*w.periods()[1]);
                auto lnsigma_real = w.ln_sigma_real(z);
                auto lnsigma_imag = w.ln_sigma_imag(z);
                auto tmp = w.sigma(z);
                if (!detail::isfinite(tmp)) {
                    continue;
                }
                auto lnsigma = std::log(tmp);
                auto err = std::abs((std::exp(lnsigma)-std::exp(complex_type(lnsigma_real,lnsigma_imag)))/std::exp(lnsigma));
                if (err > max_err) {
                    max_err = err;
                    max_g2 = g2;
                    max_g3 = g3;
                    max_z = z;
                }
                acc_err += err;
                ++counter;
            }
        }
        std::cout << "\tMax log(sigma) error: " << max_err << " @ [g2=" << max_g2 << ",g3=" << max_g3 << ",z=" << max_z << "]\n";
        std::cout << "\tAverage log(sigma) error: " << acc_err / real_type(counter) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_10)
{
    std::cout << "Testing the computation of real/imag of log(sigma)\n";
    std::cout << "==================================================\n";
    boost::mpl::for_each<real_types>(tester_10());
    std::cout << "\n\n\n";
}

struct tester_11
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        real_type g2, g3, x, P, P_comp, max_P_err = 0, acc_P_err = 0;
        size_type max_err_idx = 0;
        std::string max_err_x, max_err_P;
        for (size_type i = 0u; i < test_11_data.size(); ++i) {
            const auto &v = test_11_data[i];
            // Read the invariant values from mpmath.
            real_from_str(g2,v[0]);
            real_from_str(g3,v[1]);
            // Build the W object.
            we<real_type> w(g2,g3);
            for (decltype(v.size()) j = 2u; j < v.size(); j += 2u) {
                real_from_str(x,v[j]);
                P = w.sigma(x);
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
        std::cout << "\tMax sigma error: " << max_P_err << " @ [g2=" << test_11_data[max_err_idx][0u] << ",g3=" << test_11_data[max_err_idx][1u] << ",x=" << max_err_x << ",sigma=" << max_err_P << "]\n";
        std::cout << "\tAverage sigma error: " << acc_P_err / (real_type(test_11_data.size())*100) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_11)
{
    std::cout << "Testing the computation of real sigma\n";
    std::cout << "=====================================\n";
    boost::mpl::for_each<real_types>(tester_11());
    std::cout << "\n\n\n";
}

struct tester_12
{
    template <typename RealType>
    void operator()(const RealType &)
    {
        std::cout << "Testing type: " << typeid(RealType).name() << '\n';
        using real_type = RealType;
        std::uniform_real_distribution<double> g2dist(-10.,20.);
        std::uniform_real_distribution<double> g3dist(-10.,10.);
        real_type max_root_err = 0., max_g2 = 0., max_g3 = 0., acc_err = 0.;
        unsigned long counter = 0;
        for (int i = 0; i < 1000; ++i) {
            real_type g2(static_cast<real_type>(g2dist(rng))), g3(static_cast<real_type>(g3dist(rng)));
            we<real_type> w(g2,g3);
            auto err = std::abs((w.P(w.periods()[0].real()/real_type(2))-w.roots()[0])/w.roots()[0]);
            if (err > max_root_err) {
                max_root_err = err;
                max_g2 = g2;
                max_g3 = g3;
            }
            acc_err += err;
            ++counter;
            err = std::abs((w.P(w.periods()[1]/real_type(2))-w.roots()[2])/w.roots()[2]);
            if (err > max_root_err) {
                max_root_err = err;
                max_g2 = g2;
                max_g3 = g3;
            }
            acc_err += err;
            ++counter;
            err = std::abs((w.P((w.periods()[0]+w.periods()[1])/real_type(2))-w.roots()[1])/w.roots()[1]);
            if (err > max_root_err) {
                max_root_err = err;
                max_g2 = g2;
                max_g3 = g3;
            }
            acc_err += err;
            ++counter;
        }
        // Testing with g2 == 0.
        for (int i = 0; i < 1000; ++i) {
            real_type g2(0), g3(static_cast<real_type>(g3dist(rng)));
            we<real_type> w(g2,g3);
            auto err = std::abs((w.P(w.periods()[0].real()/real_type(2))-w.roots()[0])/w.roots()[0]);
            if (err > max_root_err) {
                max_root_err = err;
                max_g2 = g2;
                max_g3 = g3;
            }
            acc_err += err;
            ++counter;
            err = std::abs((w.P(w.periods()[1]/real_type(2))-w.roots()[2])/w.roots()[2]);
            if (err > max_root_err) {
                max_root_err = err;
                max_g2 = g2;
                max_g3 = g3;
            }
            acc_err += err;
            ++counter;
            err = std::abs((w.P((w.periods()[0]+w.periods()[1])/real_type(2))-w.roots()[1])/w.roots()[1]);
            if (err > max_root_err) {
                max_root_err = err;
                max_g2 = g2;
                max_g3 = g3;
            }
            acc_err += err;
            ++counter;
        }
        // NOTE: g3 == 0 is problematic, as one root will always be zero. Skip it for now.
        std::cout << "\tMax root match error: " << max_root_err << " @ [g2=" << max_g2 << ",g3=" << max_g3 << "]\n";
        std::cout << "\tAverage root match error: " << acc_err / real_type(counter) << '\n';
    }
};

BOOST_AUTO_TEST_CASE(test_12)
{
    std::cout << "Testing root matching\n";
    std::cout << "=====================\n";
    boost::mpl::for_each<real_types>(tester_12());
    std::cout << "\n\n\n";
}

BOOST_AUTO_TEST_CASE(test_13)
{
    we<double> w(0.62709183536928115,-0.095570490579046416);
    std::cout << w.Pinv(0.22860037855939774) << '\n';
}
