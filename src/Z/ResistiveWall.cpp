/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2017-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2017-2018: Karlsruhe Institute of Technology                 *
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
                                  , std::shared_ptr<OCLH> oclh
                                  )
  : Impedance(__calcImpedance(n,f0,f_max,L,s,xi,b),f_max,oclh)
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
