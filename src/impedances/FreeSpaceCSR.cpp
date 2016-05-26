/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Algorithms   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
 *                                                                            *
 * This file is part of Inovesa.                                              *
 * Inovesa is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Inovesa is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           *
 ******************************************************************************/

#include "impedances/FreeSpaceCSR.hpp"

vfps::FreeSpaceCSR::FreeSpaceCSR(const size_t n,
                                 const frequency_t f_rev,
                                 const frequency_t f_max)
    :
      Impedance(__calcImpedance(n,f_rev,f_max),f_max)
{
}

std::vector<vfps::impedance_t>
vfps::FreeSpaceCSR::__calcImpedance(const size_t n,
                                    const frequency_t f_rev,
                                    const frequency_t f_max)
{
    std::vector<vfps::impedance_t> rv;
    rv.reserve(n);

    // according to Eq. 6.18 in Part. Acc. Vol 57, p 35 (Murpy et al.)
    constexpr impedance_t Z0 = impedance_t(306.3,176.9);

    // frequency resolution: impedance will be sampled at multiples of delta
    const frequency_t delta = f_max/f_rev/(n-1.0);

    for (size_t i=0; i<n; i++) {
        rv.push_back(Z0*std::pow(i*delta,csrpower_t(1.0/3.0)));
    }

    return rv;
}
