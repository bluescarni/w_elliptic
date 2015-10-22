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
#include <boost/python/tuple.hpp>

#include "../src/w_elliptic.hpp"

namespace bp = boost::python;
using namespace w_elliptic;

// Wrappers for the getters.
template <typename T>
static inline bp::tuple get_invariants(const we<T> &w)
{
    return bp::make_tuple(w.invariants()[0],w.invariants()[1]);
}

template <typename T>
static inline bp::tuple get_periods(const we<T> &w)
{
    return bp::make_tuple(w.periods()[0],w.periods()[1]);
}

template <typename T>
static inline bp::tuple get_roots(const we<T> &w)
{
    return bp::make_tuple(w.roots()[0],w.roots()[1],w.roots()[2]);
}

template <typename T>
static inline T get_Delta(const we<T> &w)
{
    return w.Delta();
}

template <typename T>
static inline typename we<T>::complex_type get_q(const we<T> &w)
{
    return w.q();
}

template <typename T>
static inline bp::tuple get_etas(const we<T> &w)
{
    return bp::make_tuple(w.etas()[0],w.etas()[1]);
}

template <typename T>
static inline bp::tuple Pinv_wrapper(const we<T> &w, const typename we<T>::complex_type &c)
{
    auto res = w.Pinv(c);
    return bp::make_tuple(res[0],res[1]);
}

BOOST_PYTHON_MODULE(_core)
{
    using we_type = we<double>;
    using real_type = we_type::real_type;
    using complex_type = we_type::complex_type;
    typedef real_type (we_type::*real_1)(const real_type &) const;
    typedef complex_type (we_type::*complex_1)(const complex_type &) const;

    bp::class_<we_type> we_class("we",bp::init<const real_type &, const real_type &>());
    we_class.def(repr(bp::self));
    // First complex, then real - otherwise the Boost Python overload mechanism always
    // picks the complex part first.
    we_class.def("P",complex_1(&we_type::P));
    we_class.def("P",real_1(&we_type::P));
    we_class.def("Pprime",complex_1(&we_type::Pprime));
    we_class.def("Pprime",real_1(&we_type::Pprime));
    we_class.def("zeta",complex_1(&we_type::zeta));
    we_class.def("zeta",real_1(&we_type::zeta));
    we_class.def("sigma",complex_1(&we_type::sigma));
    we_class.def("sigma",real_1(&we_type::sigma));
    we_class.def("ln_sigma_real",&we_type::ln_sigma_real);
    we_class.def("ln_sigma_imag",&we_type::ln_sigma_imag);
    we_class.def("Pinv",Pinv_wrapper<we_type::real_type>);

    // Getters.
    we_class.add_property("invariants",get_invariants<real_type>);
    we_class.add_property("periods",get_periods<real_type>);
    we_class.add_property("roots",get_roots<real_type>);
    we_class.add_property("Delta",get_Delta<real_type>);
    we_class.add_property("q",get_q<real_type>);
    we_class.add_property("etas",get_etas<real_type>);
}
