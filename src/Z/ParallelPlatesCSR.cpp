/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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

#include "Z/ParallelPlatesCSR.hpp"

#include <boost/math/special_functions/airy.hpp>
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

vfps::ParallelPlatesCSR::ParallelPlatesCSR(const size_t n,
                                           const frequency_t f0,
                                           const frequency_t f_max,
                                           const double g)
    :
      Impedance(__calcImpedance(n,f0,f_max,g),f_max)
{
}

std::vector<vfps::impedance_t>
vfps::ParallelPlatesCSR::__calcImpedance(const size_t n,
                                         const vfps::frequency_t f0,
                                         const vfps::frequency_t f_max,
                                         const double g)
{
    std::vector<vfps::impedance_t> rv;
    rv.reserve(n);

    // frequency resolution: impedance will be sampled at multiples of delta
    const frequency_t delta = f_max/f0/(n-1.0);

    const double r_bend = physcons::c/(2*pi<double>()*f0);
    constexpr std::complex<double> j(0,1);
    rv.push_back(impedance_t(0,0));
    for (size_t i=1; i<=n/2; i++) {
        std::complex<double> Z=0;
        const double n = i*delta;
        const double m = n*std::pow(g/r_bend,3./2.);
        const uint32_t maxp = 2*m*std::pow(r_bend/g,3./2.)*f0*g/physcons::c;
        const double b = std::pow(m,-4./3.);
        std::complex<double> zinc=1;
        for (uint32_t p=1; p<=maxp;p+=2) {
            const double u = std::pow(pi<double>()*p,2)/std::pow(2,2./3.)*b;
            try {
                zinc = boost::math::airy_ai_prime(u)
                    *(boost::math::airy_ai_prime(u)
                        -j*boost::math::airy_bi_prime(u))
                    +u*boost::math::airy_ai(u)
                    *(boost::math::airy_ai(u)
                        -j*boost::math::airy_bi(u));
            } catch (...) {
                p=maxp+1;
                zinc=0;
            }
            Z += zinc;
        }
        Z *= 4.0*b*n*g/r_bend*std::pow(pi<double>(),2)*std::pow(2,1./3.)/(physcons::epsilon0*physcons::c);
        rv.push_back(impedance_t(Z));
    }
    for (size_t i=n/2+1; i<n; i++) {
        rv.push_back(impedance_t(0,0));
    }

    return rv;
}
