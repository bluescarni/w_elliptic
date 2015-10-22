.. w_elliptic documentation master file, created by
   sphinx-quickstart on Thu Oct 22 14:28:20 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to w_elliptic's documentation!
======================================

w_elliptic is a small, header-only C++11 library for the computation of `Weierstrass elliptic <https://en.wikipedia.org/wiki/Weierstrass%27s_elliptic_functions>`__
and related functions. The library has been developed with a particular focus on the application to dynamical systems,
and it offers the following features:

* it allows the definition of Weierstrassian functions in terms of real invariants :math:`g_2` and :math:`g_3`;
* the Weierstrassian functions :math:`\wp`, :math:`\wp^\prime`, :math:`\zeta`, :math:`\sigma` and :math:`\wp^{-1}`
  are supported;
* for each function (apart from :math:`\wp^{-1}`) two overloads are available, one with complex argument and one with real argument. The real
  argument overload will use only real arithmetic;
* the three standard C++ floating-point types (``float``, ``double`` and ``long double``) are supported.

Quickstart
----------

Here follows a simple usage example for the library:

.. literalinclude:: ../../tests/example.cpp
   :language: c++
   :linenos:

Installation
------------

The library does not have any dependency, apart from a C++ compiler that understands C++11. As a header-only library,
you will only need to include the ``w_elliptic.hpp`` header to start using it, no linking is needed.

Limitations
-----------

The library currently does not support neither complex invariants, nor the definition of Weierstrassian functions in terms of periods
(rather than invariants).

Contents:
=========

.. toctree::
   :maxdepth: 2

   cpp_docs.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

