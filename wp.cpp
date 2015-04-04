#include <algorithm>
#include <iomanip>
#include <boost/timer/timer.hpp>

#include "src/w_elliptic.hpp"

using namespace w_elliptic;

int main()
{
    std::cout << std::setprecision(18) << '\n';
    we<double> w(.3,.4);
    std::cout << w << '\n';
//     std::cout << w_elliptic<long double>(2.,.1) << '\n';
//     std::cout << w_elliptic<long double>(2.,-.1) << '\n';
//     std::cout << w_elliptic<long double>(.01,-.1) << '\n';
//     std::cout << w_elliptic<long double>(1E-8,3E-8) << '\n';
//     std::cout << w_elliptic<long double>(1E-8,-3E-8) << '\n';
//     std::cout << w_elliptic<long double>(1E15,-3E18) << '\n';
//     std::cout << we.P(.3) << '\n';
//     we.P(std::complex<double>(3,3));
//     return 0;
//     std::cout << we.ln_sigma_real(.1) << '\n';
//     return 0;
    {
    std::complex<double> retval = 0.;
    boost::timer::auto_cpu_timer t;
    for (int i = 0; i < 1000000; ++i) {
        retval += w.ln_sigma_real(std::complex<double>(.1+i/2000000.,.1+i/2000000.));
    }
    std::cout << retval << '\n';
    }
    std::cout << w.zeta_dup(std::complex<double>{.1,.5}) << '\n';
    auto t = w.reduce_to_fc(std::complex<double>{34.,-78.});
    std::cout << std::get<0>(t) << '\n';
    std::cout << std::get<1>(t) << '\n';
    std::cout << w.zeta(3.) << '\n';
    std::cout << "Inv:" << we<double>(1,.1).Pinv(4) << '\n';
}
