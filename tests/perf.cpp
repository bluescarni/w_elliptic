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

#include <cmath>
#include <complex>
#include <random>

#include "../src/w_elliptic.hpp"

#define BENCHPRESS_CONFIG_MAIN
#include "benchpress/src/benchpress/benchpress.hpp"

using namespace w_elliptic;

static std::mt19937 rng;

BENCHMARK("construction", [](benchpress::context* ctx) {
    ctx->stop_timer();
    std::uniform_real_distribution<double> g2dist(-10.,20.);
    std::uniform_real_distribution<double> g3dist(-10.,10.);
    for (size_t i = 0; i < ctx->num_iterations(); ++i) {
        const double g2 = g2dist(rng);
        const double g3 = g2dist(rng);
        ctx->start_timer();
        we<double> w(g2,g3);
        benchpress::escape(&w);
        ctx->stop_timer();
    }
})

BENCHMARK("real sin double", [](benchpress::context* ctx) {
    ctx->stop_timer();
    std::uniform_real_distribution<double> xdist(-10.,10.);
    for (size_t i = 0; i < ctx->num_iterations(); ++i) {
        const double x = xdist(rng);
        ctx->start_timer();
        double res = std::sin(x);
        benchpress::escape(&res);
        ctx->stop_timer();
    }
})

#define W_ELLIPTIC_BENCHMARK_REAL(Func,Type) \
BENCHMARK("real " #Func " " #Type, [](benchpress::context* ctx) { \
    ctx->stop_timer(); \
    std::uniform_real_distribution<Type> g2dist(-10.,20.); \
    std::uniform_real_distribution<Type> g3dist(-10.,10.); \
    std::uniform_real_distribution<Type> xdist(-10.,10.); \
    const unsigned nblocks = 10u; \
    const auto block_size = ctx->num_iterations() / nblocks; \
    for (auto bn = 0u; bn < nblocks; ++bn) { \
        we<Type> w(g2dist(rng),g3dist(rng)); \
        const auto block_iter = (bn == nblocks - 1u) ? (ctx->num_iterations() - bn * block_size) : block_size; \
        for (size_t i = 0; i < block_iter; ++i) { \
            const Type x = xdist(rng); \
            ctx->start_timer(); \
            auto res = w.Func(x); \
            benchpress::escape(&res); \
            ctx->stop_timer(); \
        } \
    } \
})

#define W_ELLIPTIC_BENCHMARK_COMPLEX(Func,Type) \
BENCHMARK("complex " #Func " " #Type, [](benchpress::context* ctx) { \
    ctx->stop_timer(); \
    std::uniform_real_distribution<Type> g2dist(-10.,20.); \
    std::uniform_real_distribution<Type> g3dist(-10.,10.); \
    std::uniform_real_distribution<Type> cdist(-10.,10.); \
    const unsigned nblocks = 10u; \
    const auto block_size = ctx->num_iterations() / nblocks; \
    for (auto bn = 0u; bn < nblocks; ++bn) { \
        we<Type> w(g2dist(rng),g3dist(rng)); \
        const auto block_iter = (bn == nblocks - 1u) ? (ctx->num_iterations() - bn * block_size) : block_size; \
        for (size_t i = 0; i < block_iter; ++i) { \
            const std::complex<Type> c(cdist(rng),cdist(rng)); \
            ctx->start_timer(); \
            auto res = w.Func(c); \
            benchpress::escape(&res); \
            ctx->stop_timer(); \
        } \
    } \
})

W_ELLIPTIC_BENCHMARK_REAL(P,double)
W_ELLIPTIC_BENCHMARK_REAL(Pprime,double)
W_ELLIPTIC_BENCHMARK_REAL(zeta,double)
W_ELLIPTIC_BENCHMARK_REAL(sigma,double)

BENCHMARK("complex sin double", [](benchpress::context* ctx) {
    ctx->stop_timer();
    std::uniform_real_distribution<double> cdist(-10.,10.);
    for (size_t i = 0; i < ctx->num_iterations(); ++i) {
        const std::complex<double> c(cdist(rng),cdist(rng));
        ctx->start_timer();
        auto res = std::sin(c);
        benchpress::escape(&res);
        ctx->stop_timer();
    }
})

W_ELLIPTIC_BENCHMARK_COMPLEX(P,double)
W_ELLIPTIC_BENCHMARK_COMPLEX(Pprime,double)
W_ELLIPTIC_BENCHMARK_COMPLEX(Pinv,double)
W_ELLIPTIC_BENCHMARK_COMPLEX(zeta,double)
W_ELLIPTIC_BENCHMARK_COMPLEX(sigma,double)
