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

#ifndef W_ELLIPTIC_W_ELLIPTIC_HPP
#define W_ELLIPTIC_W_ELLIPTIC_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <utility>

namespace w_elliptic
{

// Implementation details. Might be made public
// in the future.
namespace detail
{

template <typename T>
struct is_complex
{
    static const bool value = false;
};

template <typename T>
const bool is_complex<T>::value;

template <typename T>
struct is_complex<std::complex<T>>
{
    static const bool value = true;
};

template <typename T>
const bool is_complex<std::complex<T>>::value;

template <typename T, typename = void>
struct complex_t_impl
{};

template <typename T>
struct complex_t_impl<T,typename std::enable_if<std::is_floating_point<T>::value>::type>
{
    using type = std::complex<T>;
};

template <typename T, typename = void>
struct tolerance_impl
{};

template <typename T>
struct tolerance_impl<T,typename std::enable_if<std::is_floating_point<T>::value>::type>
{
    T operator()() const
    {
        return std::numeric_limits<T>::epsilon();
    }
};

template <typename T>
struct tolerance_impl<T,typename std::enable_if<is_complex<T>::value>::type>
{
    typename T::value_type operator()() const
    {
        return std::numeric_limits<typename T::value_type>::epsilon();
    }
};

template <typename T>
using complex_t = typename detail::complex_t_impl<T>::type;

template <typename T>
auto tolerance() -> decltype(detail::tolerance_impl<T>()())
{
    return detail::tolerance_impl<T>()();
}

// Number of decimal digits to print in the representation of the class.
template <typename T>
int digits10()
{
    // NOTE: print extra 2 digits - apparently this should be enough to guarantee exact
    // conversion to/from string form.
    return std::numeric_limits<T>::digits10 + 2;
}

template <typename T>
bool isfinite(const T &x)
{
    return std::isfinite(x);
}

template <typename T>
bool isfinite(const std::complex<T> &c)
{
    return isfinite(c.real()) && isfinite(c.imag());
}

}

// TODO:
// - maybe quicker calculation of error for complex: instead of abs(), just a**2+b**2;
// - testing in the fundamental cell,
// - disallow Delta = 0 for now,
// - tests with pathological values (e.g., z = 0, z close to cell border, etc.),
// - write functions that return two or three things at once (that is, P+Pprime, P+Pprime+zet),
// - comparison with Elliptic functions from Boost,
// - check the performance implications of cting floats from std::size_t in, e.g., Laurent
//   expansions,
// - proper handling of test files via cmake configure_file?
// - investigate Pprime around omega_i via Taylor expansion,
// - move setup_q stuff into sigma_setup(),
// - in series expansions, record the terms and re-accumulate them at the end in a way
//   that minimises precision loss?
template <typename T>
class we
{
    public:
        using real_type = T;
        using complex_type = detail::complex_t<real_type>;
    private:
        // TODO: checks:
        // - complex_t is defined,
        // - operations on real do not switch to other types,
        // - negation is supported,
        // - streaming is supported,
        // - construction of complex values from reals, and general interop with reals,
        // - math functions.
        static_assert(!detail::is_complex<real_type>::value,"The invariants must be reals.");
        static_assert(std::is_constructible<real_type,std::size_t>::value,"The real type must be constructible from std::size_t.");
        static_assert(std::is_constructible<real_type,int>::value,"The real type must be constructible from int.");
        static const std::size_t max_iter = 100u;
        static const real_type pi_const;
        // Coefficients for the ln sigma expansion.
        const real_type &ls_c(const std::size_t &k) const
        {
            return m_ln_sigma_c[k - 1u];
        }
        real_type &ls_c(const std::size_t &k)
        {
            return m_ln_sigma_c[k - 1u];
        }
        // Setup the roots of the Weierstrass cubic. See:
        // http://en.wikipedia.org/wiki/Cubic_function
        void calculate_roots(const real_type &a, const real_type &c, const real_type &d)
        {
            m_delta = real_type(-4) * a * c*c*c - real_type(27) * a*a * d*d;
            std::array<complex_type,3> u_list{{
                    complex_type(real_type(1)),
                    complex_type(real_type(-1),std::sqrt(real_type(3)))/real_type(2),
                    complex_type(real_type(-1),-std::sqrt(real_type(3)))/real_type(2)
            }};
            complex_type delta0 = complex_type(real_type(-3)) * complex_type(a) * complex_type(c),
                delta1 = complex_type(real_type(27)) * complex_type(a) * complex_type(a) * complex_type(d);
            // Cube root function for complex types.
            auto cbrt = [](const complex_type &c) {
                return std::pow(c,complex_type(real_type(1)/real_type(3)));
            };
            complex_type C1 = cbrt((delta1 + std::sqrt(delta1 * delta1 - complex_type(real_type(4)) * delta0 * delta0 * delta0)) / real_type(2));
            complex_type C2 = cbrt((delta1 - std::sqrt(delta1 * delta1 - complex_type(real_type(4)) * delta0 * delta0 * delta0)) / real_type(2));
            complex_type C;
            if (std::abs(C1) > std::abs(C2)) {
                C = C1;
            } else {
                C = C2;
            }
            std::transform(u_list.begin(),u_list.end(),m_roots.begin(),[&C,&delta0,&a](const complex_type &u) {
                return (real_type(-1) / (real_type(3) * a)) * (u * C + delta0 / (u * C));
            });
            if (m_delta < real_type(0)) {
                std::sort(m_roots.begin(),m_roots.end(),[](const complex_type &c1, const complex_type &c2) {
                    return c1.imag() > c2.imag();
                });
                // The second root is real.
                m_roots[1] = complex_type(m_roots[1].real(),real_type(0));
            } else {
                std::sort(m_roots.begin(),m_roots.end(),[](const complex_type &c1, const complex_type &c2) {
                    return c1.real() > c2.real();
                });
                // All the roots are real.
                std::for_each(m_roots.begin(),m_roots.end(),[](complex_type &c) {
                    c = complex_type(c.real(),real_type(0));
                });
            }
        }
        // Arithmetic-geometric mean.
        template <typename U>
        static U agm(const U &x, const U &y)
        {
            U a((x + y)/U(2)), g(std::sqrt(x*y));
            std::size_t i = 0u;
            while (true) {
                // NOTE: here it seems that the values of a/g are always around unity,
                // so we should not need to deal with relative tolerance.
                if (std::abs(a - g) <= detail::tolerance<U>()) {
                    break;
                }
                if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                U new_a = (a + g)/U(2);
                U new_g = std::sqrt(a*g);
                a = std::move(new_a);
                g = std::move(new_g);
                ++i;
            }
            return a;
        }
        // Complete elliptic integral of the first kind in terms of the modulus k.
        static complex_type K(const complex_type &k)
        {
            return pi_const / (real_type(2) * agm(real_type(1) - k,real_type(1) + k));
        }
        // See A+S 18.9. For an alternative also suitable for complex invariants, see DLMF 23.22.
        void setup_periods()
        {
            auto g2 = m_invariants[0], g3 = m_invariants[1];
            complex_type e1, e2, e3;
            e1 = m_roots[0];
            e2 = m_roots[1];
            e3 = m_roots[2];
            bool negative_g3 = false;
            if (g3 < real_type(0)) {
                // Need to change the sign of the roots
                // and invert their order too. Note that
                // delta does not change sign.
                e1 = -m_roots[2];
                e2 = -m_roots[1];
                e3 = -m_roots[0];
                g3 = -g3;
                negative_g3 = true;
            }
            complex_type om, omp, p1, p2;
            if (m_delta > real_type(0)) {
                complex_type m = (e2 - e3) / (e1 - e3);
                complex_type Km = K(std::sqrt(m));
                complex_type Kpm = K(std::sqrt(real_type(1) - m));
                om = Km / std::sqrt(e1 - e3);
                omp = complex_type(real_type(0),real_type(1)) * om * Kpm / Km;
            } else if (m_delta < real_type(0)) {
                complex_type H2(std::sqrt((e2 - e3) * (e2 - e1)));
                complex_type m(real_type(1) / real_type(2) - real_type(3) * e2 / (real_type(4) * H2));
                complex_type Km = K(std::sqrt(m));
                complex_type Kpm = K(std::sqrt(real_type(1) - m));
                complex_type om2 = Km / std::sqrt(H2);
                complex_type om2p = complex_type(real_type(0),real_type(1)) * Kpm * om2 / Km;
                om = (om2 - om2p) / real_type(2);
                omp = (om2 + om2p) / real_type(2);
            } else {
                if (g2 == real_type(0) && g3 == real_type(0)) {
                    om = complex_type(real_type(std::numeric_limits<real_type>::infinity()),real_type(0));
                    omp = complex_type(real_type(0),real_type(std::numeric_limits<real_type>::infinity()));
                } else {
                    complex_type c(e1 / real_type(2));
                    om = real_type(1) / std::sqrt(real_type(12) * c) * pi_const;
                    omp = complex_type(real_type(0),std::numeric_limits<real_type>::infinity());
                }
            }
            if (m_delta >= real_type(0)) {
                if (negative_g3) {
                    p1 = real_type(2) * omp.imag();
                    p2 = complex_type(real_type(0),real_type(2) * om.real());
                } else {
                    p1 = real_type(2) * om.real();
                    p2 = complex_type(real_type(0),real_type(2) * omp.imag());
                }
            } else {
                if (negative_g3) {
                    p1 = real_type(4) * omp.imag();
                    p2 = real_type(2) * complex_type(omp.imag(),omp.real());
                } else {
                    p1 = real_type(4) * omp.real();
                    p2 = real_type(2) * omp;
                }
            }
            m_periods[0] = p1;
            m_periods[1] = p2;
        }
        // Order the roots according to the DLMF convention. See DLMF 23.6.2 and following.
        void order_roots()
        {
            // We start by re-calculating the roots using theta functions.
            complex_type t2(m_t2[0]), t4(m_t4[0]);
            const complex_type t2_fac = real_type(2) * std::pow(m_q,real_type(1)/real_type(4));
            bool stop2 = false, stop4 = false;
            std::size_t i = 1u, counter2 = 0u, counter4 = 0u;
            while (true) {
                if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                if (stop2 && stop4) {
                    break;
                }
                if (!stop2) {
                    complex_type add2 = m_t2[i];
                    t2 += add2;
                    if (stop_check<2>(t2,add2,counter2)) {
                        stop2 = true;
                    }
                }
                if (!stop4) {
                    complex_type add4 = m_t4[i];
                    t4 += add4;
                    if (stop_check<2>(t4,add4,counter4)) {
                        stop4 = true;
                    }
                }
                ++i;
            }
            t2 *= t2_fac;
            const real_type r_fac = pi_const*pi_const/(real_type(3)*m_periods[0].real()*m_periods[0].real());
            const auto t2_4 = t2*t2*t2*t2, t4_4 = t4*t4*t4*t4;
            const complex_type e1 = r_fac * (t2_4+real_type(2)*t4_4);
            const complex_type e2 = r_fac * (t2_4-t4_4);
            const complex_type e3 = -r_fac * (real_type(2)*t2_4+t4_4);
            // Now we need to match each root calculated here to the ones we calculated (and
            // cleaned up) earlier.
            std::array<complex_type,3> roots;
            auto cmp = [](const complex_type &test, const complex_type &c1, const complex_type &c2) {
                return std::abs(test - c1) < std::abs(test - c2);
            };
            roots[0] = *std::min_element(m_roots.begin(),m_roots.end(),[&e1,cmp](const complex_type &c1, const complex_type &c2) {return cmp(e1,c1,c2);});
            roots[1] = *std::min_element(m_roots.begin(),m_roots.end(),[&e2,cmp](const complex_type &c1, const complex_type &c2) {return cmp(e2,c1,c2);});
            roots[2] = *std::min_element(m_roots.begin(),m_roots.end(),[&e3,cmp](const complex_type &c1, const complex_type &c2) {return cmp(e3,c1,c2);});
            // NOTE: maybe check that e1 is purely real, as that is used in P.
            m_roots = roots;
        }
        // Setup the q constant and related quantities.
        void setup_q()
        {
            m_q = std::pow(complex_type(-1,0),m_periods[1]/m_periods[0].real());
            // q is either pure real or pure imaginary.
            if (m_delta >= real_type(0)) {
                m_q = complex_type(m_q.real(),real_type(0));
            } else {
                m_q = complex_type(real_type(0),m_q.imag());
            }
            // Coefficients of the series expansion for ln sigma.
            real_type q2((m_q * m_q).real()), tmp(q2);
            for (std::size_t i = 1u; i < max_iter + 1u; ++i) {
                ls_c(i) = tmp/(real_type(i) * (real_type(1) - tmp));
                tmp *= q2;
            }
        }
        template <std::size_t Limit, typename U>
        static bool stop_check(const U &acc, const U &delta, std::size_t &counter)
        {
            if (!detail::isfinite(acc)) {
                return true;
            }
            if ((std::abs(acc) == real_type(0) && std::abs(delta) <= detail::tolerance<U>()) ||
                std::abs(delta/acc) <= detail::tolerance<U>())
            {
                ++counter;
                if (counter == Limit) {
                    return true;
                }
            } else {
                counter = 0u;
            }
            return false;
        }
        void setup_thetas()
        {
            const bool q_real = m_q.real() != real_type(0);
            const real_type qred = q_real ? m_q.real() : m_q.imag();
            // theta_1 and derivatives.
            for (std::size_t i = 0u; i < max_iter; ++i) {
                // TODO static check on std::size_t limits here.
                const std::size_t qpow = i*i + i;
                const real_type fac = real_type(2)*real_type(i)+real_type(1);
                const bool odd_i = i % 2u;
                real_type q(std::pow(qred,real_type(qpow)));
                // If q is pure imaginary and its power is not a multiple of 4, we need to flip the sign.
                if (!q_real && qpow % 4u) {
                    q = -q;
                }
                m_t1[i] = odd_i ? -q : q;
                m_t1p[i] = odd_i ? -q*fac : q*fac;
                m_t1ppp[i] = odd_i ? q*fac*fac*fac : -q*fac*fac*fac;
            }
            // theta_2 and derivative.
            for (std::size_t i = 0u; i < max_iter; ++i) {
                const std::size_t qpow = i*i + i;
                const real_type fac = real_type(2)*real_type(i)+real_type(1);
                real_type q(std::pow(qred,real_type(qpow)));
                if (!q_real && qpow % 4u) {
                    q = -q;
                }
                m_t2[i] = q;
                m_t2p[i] = -q*fac;
            }
            // theta_4.
            m_t4[0] = complex_type(real_type(1),real_type(0));
            for (std::size_t i = 1u; i < max_iter; ++i) {
                const std::size_t qpow = i*i;
                real_type fac = real_type(2);
                if (i % 2u) {
                    fac = -fac;
                }
                m_t4[i] = fac * std::pow(m_q,real_type(qpow));
            }
        }
        // Setup P - see DLMF 23.6.5 and round.
        void setup_P()
        {
            // We need to compute t3 * t4, which is equal to t1p/t2 according to the Jacobi identity DLMF 20.4.6.
            // This allows us to compute everything in real.
            real_type t1p(m_t1p[0]), t2(m_t2[0]);
            bool stop1 = false, stop2 = false;
            std::size_t i = 1u, counter1 = 0u, counter2 = 0u;
            while (true) {
                if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                if (stop1 && stop2) {
                    break;
                }
                if (!stop1) {
                    real_type add1 = m_t1p[i];
                    t1p += add1;
                    if (stop_check<2>(t1p,add1,counter1)) {
                        stop1 = true;
                    }
                }
                if (!stop2) {
                    real_type add2 = m_t2[i];
                    t2 += add2;
                    if (stop_check<2>(t2,add2,counter2)) {
                        stop2 = true;
                    }
                }
                ++i;
            }
            // Constant multiplicative factor for the computation of P.
            // NOTE: here t2 cannot be zero, as in order to do so tau would need to be a purely real quantity.
            m_Pfac = (pi_const * t1p/t2) / m_periods[0].real();
            m_Pfac *= m_Pfac;
        }
        // Setup Pprime.
        void setup_Pprime()
        {
            m_Pprimefac = m_Pfac * pi_const / (m_periods[0].real()/real_type(2));
        }
        void setup_zeta()
        {
            m_zetafac = pi_const/m_periods[0].real();
            m_zetaadd = m_etas[0].real()/(m_periods[0].real()/real_type(2));
        }
        // Quantities to compute sigma via theta functions. See DLMF 23.6.9.
        void setup_sigma()
        {
            // Compute the denominator of the expression of sigma via theta.
            real_type retval(m_t1p[0]);
            std::size_t i = 1u, counter = 0u;
            while (true) {
                 if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                real_type add = m_t1p[i];
                retval += add;
                if (stop_check<2>(retval,add,counter)) {
                    break;
                }
                ++i;
            }
            m_sigma_den = retval * pi_const;
        }
        // The eta constants. See DLMF 23.6.8.
        void setup_etas()
        {
            real_type t3(m_t1ppp[0]), t1(m_t1p[0]);
            bool stop3 = false, stop1 = false;
            std::size_t i = 1u, counter3 = 0u, counter1 = 0u;
            while (true) {
                if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                if (stop1 && stop3) {
                    break;
                }
                if (!stop3) {
                    real_type add3 = m_t1ppp[i];
                    t3 += add3;
                    if (stop_check<2>(t3,add3,counter3)) {
                        stop3 = true;
                    }
                }
                if (!stop1) {
                    real_type add1 = m_t1p[i];
                    t1 += add1;
                    if (stop_check<2>(t1,add1,counter1)) {
                        stop1 = true;
                    }
                }
                ++i;
            }
            m_etas[0] = complex_type(-(pi_const*pi_const)/(real_type(6)*m_periods[0].real()) * t3/t1);
            // DLMF 23.2.14.
            m_etas[1] = (m_etas[0]*m_periods[1]/real_type(2)-complex_type(0,pi_const/real_type(2)))/(m_periods[0].real()/real_type(2));
        }
        // Generic implementation of P.
        template <typename U>
        U P_impl(const U &x) const
        {
            const U arg = pi_const * x / m_periods[0].real();
            // NOTE: it seems like, in practice, GCC is smart enough to use sincos() here directly.
            const U C = std::cos(arg), S = std::sin(arg), S2 = real_type(2)*S*C, C2 = C*C-S*S;
            U Cn(C), Sn(S), tmp_s, tmp_c;
            U t2(0), t1(0);
            std::size_t i = 0u, counter2 = 0u, counter1 = 0u;
            bool stop2 = false, stop1 = false;
            while (true) {
                 if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                if (stop1 && stop2) {
                    break;
                }
                if (!stop2) {
                    U add2(m_t2[i]);
                    add2 *= Cn;
                    t2 += add2;
                    if (stop_check<2>(t2,add2,counter2)) {
                        stop2 = true;
                    }
                }
                if (!stop1) {
                    U add1(m_t1[i]);
                    add1 *= Sn;
                    t1 += add1;
                    if (stop_check<2>(t1,add1,counter1)) {
                        stop1 = true;
                    }
                }
                ++i;
                tmp_s = Sn*C2+Cn*S2;
                tmp_c = Cn*C2-Sn*S2;
                Sn = tmp_s;
                Cn = tmp_c;
            }
            return (t2*t2)/(t1*t1)*m_Pfac + m_roots[0].real();
        }
        // Generic implementation of Pprime.
        template <typename U>
        U Pprime_impl(const U &x) const
        {
            const U arg = pi_const * x / m_periods[0].real();
            const U C = std::cos(arg), S = std::sin(arg), S2 = real_type(2)*S*C, C2 = C*C-S*S;
            U Cn(C), Sn(S), tmp_s, tmp_c;
            U t1(0), t2(0), t1p(0), t2p(0);
            std::size_t i = 0u, counter1 = 0u, counter2 = 0u, counter1p = 0u, counter2p = 0u;
            bool stop1 = false, stop2 = false, stop1p = false, stop2p = false;
            while (true) {
                 if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                if (stop1 && stop2 && stop1p && stop2p) {
                    break;
                }
                if (!stop1) {
                    U add1(m_t1[i]);
                    add1 *= Sn;
                    t1 += add1;
                    if (stop_check<2>(t1,add1,counter1)) {
                        stop1 = true;
                    }
                }
                if (!stop2) {
                    U add2(m_t2[i]);
                    add2 *= Cn;
                    t2 += add2;
                    if (stop_check<2>(t2,add2,counter2)) {
                        stop2 = true;
                    }
                }
                if (!stop1p) {
                    U add1p(m_t1p[i]);
                    add1p *= Cn;
                    t1p += add1p;
                    if (stop_check<2>(t1p,add1p,counter1p)) {
                        stop1p = true;
                    }
                }
                if (!stop2p) {
                    U add2p(m_t2p[i]);
                    add2p *= Sn;
                    t2p += add2p;
                    if (stop_check<2>(t2p,add2p,counter2p)) {
                        stop2p = true;
                    }
                }
                ++i;
                tmp_s = Sn*C2+Cn*S2;
                tmp_c = Cn*C2-Sn*S2;
                Sn = tmp_s;
                Cn = tmp_c;
            }
            return m_Pprimefac*t2*(t2p*t1-t1p*t2)/(t1*t1*t1);
        }
        // Generic implementation of zeta.
        template <typename U>
        U zeta_impl(const U &x) const
        {
            const U arg = pi_const * x / m_periods[0].real();
            const U C = std::cos(arg), S = std::sin(arg), S2 = real_type(2)*S*C, C2 = C*C-S*S;
            U Cn(C), Sn(S), tmp_s, tmp_c;
            U t1(0), t1p(0);
            std::size_t i = 0u, counter1 = 0u, counter1p = 0u;
            bool stop1 = false, stop1p = false;
            while (true) {
                 if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                if (stop1 && stop1p) {
                    break;
                }
                if (!stop1) {
                    U add1(m_t1[i]);
                    add1 *= Sn;
                    t1 += add1;
                    if (stop_check<2>(t1,add1,counter1)) {
                        stop1 = true;
                    }
                }
                if (!stop1p) {
                    U add1p(m_t1p[i]);
                    add1p *= Cn;
                    t1p += add1p;
                    if (stop_check<2>(t1p,add1p,counter1p)) {
                        stop1p = true;
                    }
                }
                ++i;
                tmp_s = Sn*C2+Cn*S2;
                tmp_c = Cn*C2-Sn*S2;
                Sn = tmp_s;
                Cn = tmp_c;
            }
            return m_zetaadd*x+m_zetafac*t1p/t1;
        }
        // Generic sigma implementation.
        template <typename U>
        U sigma_impl(const U &x) const
        {
            const U arg = pi_const * x / m_periods[0].real();
            const U C = std::cos(arg), S = std::sin(arg), S2 = real_type(2)*S*C, C2 = C*C-S*S;
            U Cn(C), Sn(S), tmp_s, tmp_c;
            U retval(0);
            std::size_t i = 0u, counter = 0u;
            while (true) {
                 if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                U add(m_t1[i]);
                add *= Sn;
                retval += add;
                if (stop_check<2>(retval,add,counter)) {
                    break;
                }
                ++i;
                tmp_s = Sn*C2+Cn*S2;
                tmp_c = Cn*C2-Sn*S2;
                Sn = tmp_s;
                Cn = tmp_c;
            }
            // The rest.
            retval /= m_sigma_den;
            retval *= m_periods[0].real();
            retval *= std::exp(m_etas[0].real()*x*x/m_periods[0].real());
            return retval;
        }
        // Continuous implementations of the real/imag parts of ln sigma.
        real_type ln_sigma_imag_cont(const complex_type &c) const
        {
            const real_type om1 = m_periods[0].real()/real_type(2);
            // Reduce to the fundamental real period.
            real_type N(std::floor(c.real()/m_periods[0].real()));
            complex_type cred_F(c.real() - N*m_periods[0].real(),c.imag());
            // Perform the expansion.
            complex_type arg = pi_const * cred_F / m_periods[0].real();
            const real_type a = arg.real(), b = arg.imag();
            real_type retval = m_etas[0].real() * real_type(2) * cred_F.real()* cred_F.imag() / m_periods[0].real()
                + std::log(std::sin(pi_const*cred_F/m_periods[0].real())).imag();
            std::size_t i = 1u, miter = max_iter + 1u, counter = 0u;
            // NOTE: here the sincos function might be useful.
            real_type C = std::cos(a), S = std::sin(a), Cn(1), Sn(0), Ch = std::cosh(b), Sh = std::sinh(b), Chn(1), Shn(0);
            real_type tmp_s, tmp_c, tmp_sh, tmp_ch, mul;
            while (true) {
                if (i == miter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                // These are just the addition formulae for trig and hyperbolic functions.
                tmp_s = Sn*C+Cn*S;
                tmp_c = Cn*C-Sn*S;
                tmp_sh = Shn*Ch+Chn*Sh;
                tmp_ch = Chn*Ch+Shn*Sh;
                Sn = tmp_s;
                Cn = tmp_c;
                Shn = tmp_sh;
                Chn = tmp_ch;
                mul = real_type(8)*Sn*Cn*Shn*Chn;
                real_type add = ls_c(i) * mul;
                retval += add;
                if (stop_check<2>(retval,add,counter)) {
                    break;
                }
                ++i;
            }
            // Add the homogeneity relations.
            retval += real_type(2)*N*m_etas[0].real()*cred_F.imag() - N*pi_const;
            return retval;
        }
        real_type ln_sigma_real_cont(const complex_type &c) const
        {
            const real_type om1 = m_periods[0].real()/real_type(2);
            // Reduce to the fundamental real period.
            real_type N(std::floor(c.real()/m_periods[0].real()));
            complex_type cred_F(c.real() - N*m_periods[0].real(),c.imag());
            // Perform the expansion.
            complex_type arg = pi_const * cred_F / m_periods[0].real();
            const real_type a = arg.real(), b = arg.imag();
            real_type retval = std::log(m_periods[0].real()/pi_const) + m_etas[0].real() * (cred_F.real()*cred_F.real()-cred_F.imag()*cred_F.imag()) / m_periods[0].real()
                + std::log(std::sin(pi_const*cred_F/m_periods[0].real())).real();
            std::size_t i = 1u, miter = max_iter + 1u, counter = 0u;
            // NOTE: here the sincos function might be useful.
            real_type C = std::cos(a), S = std::sin(a), Cn(1), Sn(0), Ch = std::cosh(b), Sh = std::sinh(b), Chn(1), Shn(0);
            real_type tmp_s, tmp_c, tmp_sh, tmp_ch, mul;
            while (true) {
                if (i == miter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                tmp_s = Sn*C+Cn*S;
                tmp_c = Cn*C-Sn*S;
                tmp_sh = Shn*Ch+Chn*Sh;
                tmp_ch = Chn*Ch+Shn*Sh;
                Sn = tmp_s;
                Cn = tmp_c;
                Shn = tmp_sh;
                Chn = tmp_ch;
                mul = Sn*Sn*Chn*Chn-Cn*Cn*Shn*Shn;
                mul *= real_type(4);
                real_type add = ls_c(i) * mul;
                retval += add;
                if (stop_check<2>(retval,add,counter)) {
                    break;
                }
                ++i;
            }
            // Add the homogeneity relations.
            retval += real_type(2)*N*m_etas[0].real()*(cred_F.real() + N*om1);
            return retval;
        }
        // Reduction to fundamental cell.
        std::tuple<real_type,real_type> reduce_to_fc(const complex_type &c) const
        {
            real_type re(c.real()), im(c.imag());
            real_type p1r(m_periods[0].real()), p1i(m_periods[0].imag()), p2r(m_periods[1].real()), p2i(m_periods[1].imag());
            real_type det(p1r*p2i-p2r*p1i);
            real_type alpha((p2i*re-p2r*im)/det), beta((-p1i*re+p1r*im)/det);
            return std::make_tuple(std::move(alpha),std::move(beta));
        }
    public:
        // TODO:
        // - finiteness checks,
        // - handle special cases -> singularities in the cubic roots computations, infinite periods,
        //   stuff like that.
        explicit we(const real_type &g2, const real_type &g3):m_invariants{{g2,g3}}
        {
            // Calculation of the roots. Note that after this the roots
            // are stored in an unspecified order.
            calculate_roots(real_type(4),-g2,-g3);
            // Computation of the periods.
            setup_periods();
            // Setup q and related.
            setup_q();
            // Setup of the needed theta functions expansions.
            setup_thetas();
            // After setting up q and the theta expansion, we can order the roots according to
            // the DLMF convention.
            order_roots();
            // Calculate the eta constants.
            setup_etas();
            // Setup the quantities for the calculation of P.
            setup_P();
            // Setup the quantities for the calculation of Pprime.
            setup_Pprime();
            // Setup the quantities for the calculation of zeta.
            setup_zeta();
            // Setup the quantities for the computation of sigma.
            setup_sigma();
        }
        const std::array<complex_type,3> &roots() const
        {
            return m_roots;
        }
        const std::array<real_type,2> &invariants() const
        {
            return m_invariants;
        }
        const std::array<complex_type,2> &periods() const
        {
            return m_periods;
        }
        const std::array<complex_type,2> &etas() const
        {
            return m_etas;
        }
        const real_type &Delta() const
        {
            return m_delta;
        }
        const complex_type &q() const
        {
            return m_q;
        }
        friend std::ostream &operator<<(std::ostream &os, const we &w)
        {
            std::ostringstream oss;
            oss << std::setprecision(detail::digits10<real_type>());
            oss << "Invariants: [" << w.m_invariants[0] << ',' << w.m_invariants[1] << "]\n";
            oss << "Delta: " << w.m_delta << '\n';
            oss << "Roots: [" << w.m_roots[0] << ',' << w.m_roots[1] << ',' << w.m_roots[2] << "]\n";
            oss << "Periods: [" << w.m_periods[0] << ',' << w.m_periods[1] << "]\n";
            oss << "etas: [" << w.m_etas[0] << ',' << w.m_etas[1] << "]\n";
            oss << "q: " << w.m_q;
            os << oss.str();
            return os;
        }
        real_type P(const real_type &x) const
        {
            return P_impl(x);
        }
        complex_type P(const complex_type &c) const
        {
            auto ab = reduce_to_fc(c);
            real_type alpha = std::get<0>(ab) - std::floor(std::get<0>(ab)),
                beta = std::get<1>(ab) - std::floor(std::get<1>(ab));
            complex_type cred = m_periods[0].real() * alpha + m_periods[1] * beta;
            return P_impl(cred);
        }
        real_type Pprime(const real_type &x) const
        {
            return Pprime_impl(x);
        }
        complex_type Pprime(const complex_type &c) const
        {
            auto ab = reduce_to_fc(c);
            real_type alpha = std::get<0>(ab) - std::floor(std::get<0>(ab)),
                beta = std::get<1>(ab) - std::floor(std::get<1>(ab));
            complex_type cred = m_periods[0].real() * alpha + m_periods[1] * beta;
            return Pprime_impl(cred);
        }
        real_type zeta(const real_type &x) const
        {
            return zeta_impl(x);
        }
        complex_type zeta(const complex_type &c) const
        {
            // Reduction to the fundamental cell.
            auto ab = reduce_to_fc(c);
            real_type N = std::floor(std::get<0>(ab)), M = std::floor(std::get<1>(ab));
            real_type alpha = std::get<0>(ab) - N, beta = std::get<1>(ab) - M;
            complex_type cred(m_periods[0].real() * alpha + m_periods[1] * beta);
            return zeta_impl(cred) + real_type(2)*N*m_etas[0] + real_type(2)*M*m_etas[1];
        }
        real_type sigma(const real_type &x) const
        {
            return sigma_impl(x);
        }
        complex_type sigma(const complex_type &c) const
        {
            return sigma_impl(c);
        }
        real_type ln_sigma_real(const complex_type &c) const
        {
            if (c.imag() >= real_type(0) && c.imag() <= m_periods[1].imag()/real_type(2)) {
                return ln_sigma_real_cont(c);
            }
            return std::log(sigma(c)).real();
        }
        real_type ln_sigma_imag(const complex_type &c) const
        {
            if (c.imag() >= real_type(0) && c.imag() <= m_periods[1].imag()/real_type(2)) {
                return ln_sigma_imag_cont(c);
            }
            return std::log(sigma(c)).imag();
        }
        // Inversion of DLMF 23.6.21.
        complex_type Pinv(const complex_type &c) const
        {
            complex_type e1, e2, e3;
            e1 = m_roots[0];
            e2 = m_roots[1];
            e3 = m_roots[2];
            // Check if the computation will have a singularity.
            bool singular = false;
            if (e3 == real_type(0)) {
                if (std::abs(c - e3) <= detail::tolerance<real_type>()) {
                    singular = true;
                }
            } else {
                if (std::abs((c - e3)/e3) <= detail::tolerance<real_type>()) {
                    singular = true;
                }
            }
            complex_type k = std::sqrt((e2-e3)/(e1-e3)), phi(singular ? complex_type(real_type(0)) : std::asin(std::sqrt((e1-e3)/(c-e3)))),
                tmp(real_type(1));
            std::size_t i = 0u;
            while (true) {
                if (i == max_iter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                if (std::abs(k - real_type(1)) <= detail::tolerance<complex_type>()) {
                    break;
                }
                if (!singular) {
                    phi = (std::asin(k*std::sin(phi)) + phi)/real_type(2);
                }
                tmp *= real_type(2) / (real_type(1) + k);
                k = real_type(2) * std::sqrt(k)  /(real_type(1) + k);
                ++i;
            }
            complex_type ell;
            if (singular) {
                // NOTE: strictly speaking, the sign of the factor here depends
                // on the sign of the infinity generated above. But we are going
                // to reduce to the fpp anyway, so the sign does not really matter.
                ell = complex_type(0,pi_const/real_type(2)) * tmp;
            } else {
                ell = std::log(real_type(1)/std::cos(phi)+std::tan(phi)) * tmp;
            }
            complex_type retval = ell / std::sqrt(e1-e3);
            // Reduction to the fundamental cell.
            auto ab = reduce_to_fc(retval);
            real_type alpha = std::get<0>(ab) - std::floor(std::get<0>(ab)),
                beta = std::get<1>(ab) - std::floor(std::get<1>(ab));
            retval = alpha * m_periods[0].real() + beta * m_periods[1];
            // Pick the value with the smallest imaginary part.
            auto alt_retval = -retval + m_periods[0].real() + m_periods[1];
            if (alt_retval.imag() < retval.imag()) {
                retval = alt_retval;
            }
            return retval;
        }
    private:
        std::array<real_type,2>                     m_invariants;
        real_type                                   m_delta;
        std::array<complex_type,3>                  m_roots;
        std::array<complex_type,2>                  m_periods;
        std::array<complex_type,2>                  m_etas;
        complex_type                                m_q;
        std::array<real_type,max_iter>              m_ln_sigma_c;
        std::array<real_type,max_iter>              m_t1;
        std::array<real_type,max_iter>              m_t1p;
        std::array<real_type,max_iter>              m_t1ppp;
        std::array<real_type,max_iter>              m_t2;
        std::array<real_type,max_iter>              m_t2p;
        std::array<complex_type,max_iter>           m_t4;
        // Constants for the computations of the various functions.
        real_type                                   m_Pfac;
        real_type                                   m_Pprimefac;
        real_type                                   m_sigma_den;
        real_type                                   m_zetafac;
        real_type                                   m_zetaadd;
};

template <typename T>
const typename we<T>::real_type we<T>::pi_const = std::acos(typename we<T>::real_type(-1));

template <typename T>
const std::size_t we<T>::max_iter;

}

#endif
