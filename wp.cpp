#include <algorithm>
#include <array>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>

#include <iomanip>
#include <boost/timer/timer.hpp>

namespace weierstrass_elliptic
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

namespace detail
{

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

}

template <typename T>
using complex_t = typename detail::complex_t_impl<T>::type;

template <typename T>
auto tolerance() -> decltype(detail::tolerance_impl<T>()())
{
    return detail::tolerance_impl<T>()();
}

template <typename T>
class w_elliptic
{
    public:
        using real_type = T;
        using complex_type = complex_t<real_type>;
    private:
        // TODO: checks:
        // - complex_t is defined,
        // - operations on real do not switch to other types,
        // - negation is supported,
        // - streaming is supported,
        // - construction of complex values from reals, and general interop with reals,
        // - math functions.
        static_assert(!is_complex<real_type>::value,"The invariants must be reals.");
        static_assert(std::is_constructible<real_type,std::size_t>::value,"The real type must be constructible from std::size_t.");
        static_assert(std::is_constructible<real_type,int>::value,"The real type must be constructible from int.");
        static const std::size_t max_iter = 100u;
        // Utilities to access the coefficients ck of the Laurent series expansions
        // for the Weierstrassian functions.
        const real_type &ck(const std::size_t &k) const
        {
            return m_ck_laurent[k - 2u];
        }
        real_type &ck(const std::size_t &k)
        {
            return m_ck_laurent[k - 2u];
        }
        // Setup the roots of the Weierstrass cubic. See:
        // http://en.wikipedia.org/wiki/Cubic_function
        void setup_roots(const real_type &a, const real_type &c, const real_type &d)
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
                if (std::abs(a - g) <= tolerance<U>()) {
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
        // Setup of the Laurent series coefficients: A+S 18.5.
        void setup_laurent()
        {
            // Construct the coefficients of the Laurent expansions for P, P', etc.
            // This will be c2.
            ck(2u) = m_invariants[0]/real_type(20);
            // c3.
            ck(3u) = m_invariants[1]/real_type(28);
            // The rest: A+S 18.5.3.
            for (std::size_t k = 4u; k < max_iter + 2u; ++k) {
                // Accumulator.
                real_type acc(ck(2u) * ck(k - 2u));
                for (std::size_t m = 3u; m <= k - 2u; ++m) {
                    acc += ck(m) * ck(k - m);
                }
                acc *= real_type(3) / ((real_type(2)*real_type(k) + real_type(1)) * (real_type(k) - real_type(3)));
                ck(k) = acc;
            }
        }
        // Complete elliptic integral of the first kind in terms of the modulus k.
        static complex_type K(const complex_type &k)
        {
            return boost::math::constants::pi<real_type>() / (real_type(2) * agm(real_type(1) - k,real_type(1) + k));
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
                    om = real_type(1) / std::sqrt(real_type(12) * c) * boost::math::constants::pi<real_type>();
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
            // Store also the radius of convergence of the Laurent series of P, Pprime and zeta.
            // See note at the end of A+S 18.5.
            std::array<real_type,4> candidates = {{std::abs(p1),std::abs(p2),std::abs(p1+p2),std::abs(p1-p2)}};
            m_conv_radius = *std::min_element(candidates.begin(),candidates.end());
// std::cout << "m_conv_radius:" << m_conv_radius << '\n';
        }
        // Duplication formula for P.
        template <typename U, typename V>
        static U duplicate_P(const U &Px, const V &g2, const V &g2_2, const V &g3)
        {
            U P2(Px*Px), P3 = (P2*Px), num(V(6)*P2 - g2_2);
            num *= num;
            return V(-2)*Px + num / (V(4)*(V(4)*P3-g2*Px-g3));
        }
        // Duplication formula for Pprime.
        template <typename U, typename V>
        static std::tuple<U,U> duplicate_Pprime(const U &Ppx, const U &Px, const V &g2, const V &g2_2, const V &g3)
        {
            U Pp2(Ppx*Ppx), Pp3(Pp2*Ppx), Pp4(Pp2*Pp2), Ppp(V(6)*Px*Px-g2_2), Ppp3(Ppp*Ppp*Ppp);
            return std::make_tuple((V(-4)*Pp4+V(12)*Px*Pp2*Ppp-Ppp3)/(V(4)*Pp3),duplicate_P(Px,g2,g2_2,g3));
        }
        // Duplication formula for zeta.
        template <typename U, typename V>
        static std::tuple<U,U,U> duplicate_zeta(const U &zetax, const U &Ppx, const U &Px, const V &g2, const V &g2_2, const V &g3)
        {
            U Ppp(V(6)*Px*Px-g2_2);
            auto dup_Pp = duplicate_Pprime(Ppx,Px,g2,g2_2,g3);
            return std::make_tuple(U(2)*zetax+Ppp/(U(2)*Ppx),std::get<0>(dup_Pp),std::get<1>(dup_Pp));
        }
    public:
        // TODO:
        // - finiteness checks,
        // - handle special cases -> singularities in the cubic roots computations, infinite periods,
        //   stuff like that.
        explicit w_elliptic(const real_type &g2, const real_type &g3):m_invariants{{g2,g3}}
        {
            // Setup Laurent series coefficients.
            setup_laurent();
            // Calculation of the roots.
            setup_roots(real_type(4),-g2,-g3);
            // Computation of the periods.
            setup_periods();
        }
        const std::array<real_type,2> &invariants() const
        {
            return m_invariants;
        }
        const std::array<complex_type,2> &periods() const
        {
            return m_periods;
        }
        friend std::ostream &operator<<(std::ostream &os, const w_elliptic &w)
        {
            os << "Invariants: [" << w.m_invariants[0] << ',' << w.m_invariants[1] << "]\n";
            os << "Delta: " << w.m_delta << '\n';
            os << "Roots: [" << w.m_roots[0] << ',' << w.m_roots[1] << ',' << w.m_roots[2] << "]\n";
            os << "Periods: [" << w.m_periods[0] << ',' << w.m_periods[1] << ']';
            return os;
        }
        real_type P(const real_type &x) const
        {
            const real_type g2 = m_invariants[0], g3 = m_invariants[1], g2_2 = g2/real_type(2);
            real_type xred(x);
            // P is an even function.
            if (xred < real_type(0)) {
                xred = -xred;
            }
            // Reduction of x to the fundamental cell.
            xred = std::fmod(xred,m_periods[0].real());
            // Further reduction.
            if (xred > m_periods[0].real() / real_type(2)) {
                xred = m_periods[0].real() - xred;
            }
            // Now we need to reduce xred to the radius of convergence of the Laurent series.
            std::size_t n = 0u;
            while (xred >= m_conv_radius / real_type(8)) {
                xred /= real_type(2);
                ++n;
            }
            real_type retval(P_laurent(xred));
            for (std::size_t i = 0u; i < n; ++i) {
                retval = duplicate_P(retval,g2,g2_2,g3);
            }
            return retval;
        }
        real_type Pprime(const real_type &x) const
        {
            const real_type g2 = m_invariants[0], g3 = m_invariants[1], g2_2 = g2/real_type(2);
            real_type xred(x);
            // Pprime is an odd function.
            bool negate = false;
            if (xred < real_type(0)) {
                xred = -xred;
                negate = !negate;
            }
std::cout << "xred:" << xred << '\n';
            // Reduction of x to the fundamental cell.
            xred = std::fmod(xred,m_periods[0].real());
            // Further reduction.
            if (xred > m_periods[0].real() / real_type(2)) {
                xred = m_periods[0].real() - xred;
                negate = !negate;
            }
            // Now we need to reduce xred to the radius of convergence of the Laurent series.
            std::size_t n = 0u;
            while (xred >= m_conv_radius / real_type(8)) {
                xred /= real_type(2);
                ++n;
            }
            auto retval = std::make_tuple(Pprime_laurent(xred),P_laurent(xred));
            for (std::size_t i = 0u; i < n; ++i) {
                retval = duplicate_Pprime(std::get<0>(retval),std::get<1>(retval),g2,g2_2,g3);
            }
            if (negate) {
                return -std::get<0>(retval);
            }
            return std::get<0>(retval);
        }
        real_type zeta_dup(const real_type &x) const
        {
            const real_type g2 = m_invariants[0], g3 = m_invariants[1], g2_2 = g2/real_type(2);
            real_type xred(x);
            // We need to reduce x to the radius of convergence of the Laurent series.
            std::size_t n = 0u;
            while (xred >= m_conv_radius / real_type(8)) {
                xred /= real_type(2);
                ++n;
            }
            auto retval = std::make_tuple(zeta_laurent(xred),Pprime_laurent(xred),P_laurent(xred));
            for (std::size_t i = 0u; i < n; ++i) {
                retval = duplicate_zeta(std::get<0>(retval),std::get<1>(retval),std::get<2>(retval),g2,g2_2,g3);
            }
            return std::get<0>(retval);
        }
        template <typename U>
        U P_laurent(const U &z) const
        {
            U z2(z*z), tmp(z2);
            U retval = U(1) / z2;
            std::size_t i = 2u, miter = max_iter + 2u;
            while (true) {
                if (i == miter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                U add(ck(i) * tmp);
                retval += add;
                if (std::abs(add/retval) <= tolerance<U>()) {
                    break;
                }
                tmp *= z2;
                ++i;
            }
            return retval;
        }
        template <typename U>
        U Pprime_laurent(const U &z) const
        {
            U z2(z*z), tmp(z);
            U retval = U(-2) / (z2*z);
            std::size_t i = 2u, miter = max_iter + 2u;
            while (true) {
                if (i == miter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                U add((U(2)*U(i) - U(2))*ck(i) * tmp);
                retval += add;
                if (std::abs(add/retval) <= tolerance<U>()) {
                    break;
                }
                tmp *= z2;
                ++i;
            }
            return retval;
        }
        template <typename U>
        U zeta_laurent(const U &z) const
        {
            U z2(z*z), tmp(z*z2);
            U retval = U(1) / (z);
            std::size_t i = 2u, miter = max_iter + 2u;
            while (true) {
                if (i == miter) {
                    std::cout << "WARNING max_iter reached\n";
                    break;
                }
                U sub(ck(i) * tmp/(U(2)*U(i) - U(1)));
                retval -= sub;
                if (std::abs(sub/retval) <= tolerance<U>()) {
                    break;
                }
                tmp *= z2;
                ++i;
            }
            return retval;
        }
    private:
        std::array<real_type,2>         m_invariants;
        real_type                       m_delta;
        std::array<complex_type,3>      m_roots;
        std::array<real_type,max_iter>  m_ck_laurent;
        std::array<complex_type,2>      m_periods;
        real_type                       m_conv_radius;
};

}

using namespace weierstrass_elliptic;

int main()
{
    std::cout << std::setprecision(16) << '\n';
    w_elliptic<double> we(.3,.4);
//     std::cout << we << '\n';
//     std::cout << w_elliptic<long double>(2.,.1) << '\n';
//     std::cout << w_elliptic<long double>(2.,-.1) << '\n';
//     std::cout << w_elliptic<long double>(.01,-.1) << '\n';
//     std::cout << w_elliptic<long double>(1E-8,3E-8) << '\n';
//     std::cout << w_elliptic<long double>(1E-8,-3E-8) << '\n';
//     std::cout << w_elliptic<long double>(1E15,-3E18) << '\n';
//     std::cout << we.P(.3) << '\n';
//     double retval = 0.;
//     boost::timer::auto_cpu_timer t;
//     for (int i = 0; i < 1000; ++i) {
//         retval += we.Pprime(2.+i/1000.);
//     }
//     std::cout << retval << '\n';
    std::cout << we.zeta_dup(2.) << '\n';

}
