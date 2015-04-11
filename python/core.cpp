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

// NOTE: the order of inclusion in the first two items here is forced by these two issues:
// http://mail.python.org/pipermail/python-list/2004-March/907592.html
// http://mail.python.org/pipermail/new-bugs-announce/2011-March/010395.html
#if defined(_WIN32)
#include <cmath>
#include <Python.h>
#else
#include <Python.h>
#include <cmath>
#endif

#if PY_MAJOR_VERSION < 2 || (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 6)
    #error Minimum supported Python version is 2.6.
#endif

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>

#include "../src/w_elliptic.hpp"

namespace bp = boost::python;
using namespace w_elliptic;

BOOST_PYTHON_MODULE(_core)
{
    using we_type = we<double>;
    using real_type = we_type::real_type;
    using complex_type = we_type::complex_type;
    typedef real_type (we_type::*real_1)(const real_type &) const;
    typedef complex_type (we_type::*complex_1)(const complex_type &) const;

    bp::class_<we_type> we_class("we",bp::init<const double &, const double &>());
    we_class.def(repr(bp::self));
    // First complex, then real - otherwise the Boost Python overload mechanism always
    // picks the complex part first.
    we_class.def("P",complex_1(&we_type::P));
    we_class.def("P",real_1(&we_type::P));
    we_class.def("Pprime",complex_1(&we_type::Pprime));
    we_class.def("Pprime",real_1(&we_type::Pprime));
    we_class.def("zeta",complex_1(&we_type::zeta));
    we_class.def("zeta",real_1(&we_type::zeta));
    we_class.def("sigma",&we_type::sigma);
    we_class.def("ln_sigma",&we_type::ln_sigma);
    we_class.def("ln_sigma_real_cont",&we_type::ln_sigma_real_cont);
    we_class.def("ln_sigma_imag_cont",&we_type::ln_sigma_imag_cont);
}
