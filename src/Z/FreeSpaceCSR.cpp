// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include "Z/FreeSpaceCSR.hpp"

vfps::FreeSpaceCSR::FreeSpaceCSR( const size_t n
                                , const frequency_t f_rev
                                , const frequency_t f_max
                                , oclhptr_t oclh
                                )
  : Impedance( __calcImpedance(n,f_rev,f_max),f_max, oclh)
{
}

std::vector<vfps::impedance_t>
vfps::FreeSpaceCSR::__calcImpedance(const size_t n,
                                    const frequency_t f_rev,
                                    const frequency_t f_max)
{
    std::vector<vfps::impedance_t> rv;
    rv.reserve(n);

    /* Zeros contribution to the impedance
     * according to Eq. 6.18 in Part. Acc. Vol 57, p 35 (Murpy et al.)
     * (Note that this is not the free space impedance.)
     */
    constexpr impedance_t Z0 = impedance_t(306.3,176.9);

    // frequency resolution: impedance will be sampled at multiples of delta
    const frequency_t delta = f_max/f_rev/(n-1.0);

    for (size_t i=0; i<=n/2; i++) {
        rv.push_back(Z0*std::pow(i*delta,csrpower_t(1.0/3.0)));
    }
    for (size_t i=n/2+1; i<n; i++) {
        rv.push_back(impedance_t(0,0));
    }

    return rv;
}
