// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#include "Z/ParallelPlatesCSR.hpp"

#include <boost/math/special_functions/airy.hpp>
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;
using boost::math::constants::pi_sqr;

vfps::ParallelPlatesCSR::ParallelPlatesCSR( const size_t nfreqs
                                          , const frequency_t f0
                                          , const frequency_t f_max
                                          , const double g
                                          , oclhptr_t oclh
                                          )
    :
      Impedance(__calcImpedance(nfreqs,f0,f_max,g),f_max, oclh )
{
}

std::vector<vfps::impedance_t>
vfps::ParallelPlatesCSR::__calcImpedance(const size_t nfreqs,
                                         const vfps::frequency_t f0,
                                         const vfps::frequency_t f_max,
                                         const double g)
{
    std::vector<vfps::impedance_t> rv(nfreqs,0);

    // frequency resolution: impedance will be sampled at multiples of delta
    const frequency_t delta = f_max/f0/(nfreqs-1.0);

    const double r_bend = physcons::c/(2*pi<double>()*f0);
    constexpr std::complex<double> j(0,1);
    for (size_t i=1; i<=nfreqs/2; i++) {
        std::complex<double> Z=0;
        const double n = i*delta;
        const double m = n*std::pow(g/r_bend,3./2.);
        const uint32_t maxp = 2*m*std::pow(r_bend/g,3./2.)*f0*g/physcons::c;
        const double b = std::pow(m,-4./3.);
        std::complex<double> zinc;
        for (uint32_t p=1; p<=maxp;p+=2) {
            const double u = pi_sqr<double>()*p*p/std::pow(2,2./3.)*b;
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
        Z *= 4.0*b*n*g/r_bend*pi_sqr<double>()*std::pow(2,1./3.)*physcons::Z0;
        rv[i] = Z;
    }

    return rv;
}
