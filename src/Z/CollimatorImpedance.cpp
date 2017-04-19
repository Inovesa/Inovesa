/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2017: Patrik Schönfeldt                                      *
 * Copyright (c) 2017: Karlsruhe Institute of Technology                      *
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

#include "Z/CollimatorImpedance.hpp"

vfps::CollimatorImpedance::CollimatorImpedance(const size_t n,
                                   const frequency_t f_max,
                                   const double outer,
                                   const double inner)
    :
      Impedance(__calcImpedance(n,outer,inner),f_max)
{
}

std::vector<vfps::impedance_t>
vfps::CollimatorImpedance::__calcImpedance(const size_t n,
                                     const double outer,
                                     const double inner)
{
    std::vector<vfps::impedance_t> rv;
    rv.reserve(n);

    const impedance_t Z(Z0/M_PI*std::log(outer/inner),0);

    for (size_t i=0; i<=n/2; i++) {
        rv.push_back(Z);
    }
    for (size_t i=n/2+1; i<n; i++) {
        rv.push_back(impedance_t(0,0));
    }

    return rv;
}
