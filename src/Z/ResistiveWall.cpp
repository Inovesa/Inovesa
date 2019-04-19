// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include "Z/ResistiveWall.hpp"
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

vfps::ResistiveWall::ResistiveWall( const size_t n
                                  , const frequency_t f0
                                  , const frequency_t f_max
                                  , const double L
                                  , const double s
                                  , const double xi
                                  , const double b
                                  , oclhptr_t oclh
                                  )
  : Impedance( __calcImpedance(n,f0,f_max,L,s,xi,b),f_max, oclh)
{
}

std::vector<vfps::impedance_t>
vfps::ResistiveWall::__calcImpedance(const size_t n,
                                     const frequency_t f0,
                                     const frequency_t f_max,
                                     const double L,
                                     const double s,
                                     const double xi,
                                     const double b)
{
    std::vector<vfps::impedance_t> rv;
    rv.reserve(n);

    const double mu_r = (1+xi);

    /* First contribution to the impedance (at n=1).
     * Note that we do not care for n<0, so we can drop the sgn(n).
     */
    const impedance_t Z1 =
            static_cast<frequency_t>(
                std::sqrt(Z0*mu_r*f0/s/pi<double>()/physcons::c)*L/2/b
            ) * impedance_t(1,-1);

    // frequency resolution: impedance will be sampled at multiples of delta
    const frequency_t delta = f_max/f0/(n-1.0);

    for (size_t i=0; i<=n/2; i++) {
        rv.push_back(Z1*std::sqrt(i*delta));
    }
    for (size_t i=n/2+1; i<n; i++) {
        rv.push_back(impedance_t(0,0));
    }

    return rv;
}
