// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#include "Z/ConstImpedance.hpp"

vfps::ConstImpedance::ConstImpedance( const size_t n
                                    , const frequency_t f_max
                                    , const vfps::impedance_t Z
                                    , oclhptr_t oclh
                                    )
  : Impedance( __calcImpedance(n,Z),f_max, oclh)
{
}

std::vector<vfps::impedance_t>
vfps::ConstImpedance::__calcImpedance(const size_t n,
                                      const vfps::impedance_t Z)
{
    std::vector<vfps::impedance_t> rv;
    rv.reserve(n);
    rv.resize(n/2,Z);
    rv.resize(n,0);

    return rv;
}
