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
 *along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.            *
 ******************************************************************************/

#include "Impedance.hpp"

#include <boost/math/special_functions/airy.hpp>

vfps::Impedance::Impedance(vfps::Impedance::ImpedanceModel model, size_t n,
                           const double f_rev, const double f_max) :
    _nmax(n),
    _axis(Ruler<frequency_t>(_nmax,0,f_max,1))
{
    _data.reserve(_nmax);

    // according to Eq. 6.18 in Part. Acc. Vol 57, p 35 (Murpy et al.)
    constexpr impedance_t Z0 = impedance_t(306.3,176.9);

    // frequency resolution: impedance will be sampled at multiples of delta
    const frequency_t delta = f_max/f_rev/(_nmax-1.0);

    switch (model) {
    case ImpedanceModel::ConstantOne:
        for (size_t i=0; i<_nmax; i++) {
            _data.push_back(impedance_t(1,0));
        }
        break;
    case ImpedanceModel::FreeSpaceCSR:
        for (size_t i=0; i<_nmax; i++) {
            _data.push_back(Z0*std::pow(i*delta,csrpower_t(1.0/3.0)));
        }
        break;
    case ImpedanceModel::ParallelPlates:
        uint32_t maxp = 3;
        csrpower_t h_r = 0.7;
        constexpr impedance_t j(0,1);
        for (size_t i=0; i<_nmax; i++) {
            impedance_t Z=0;
            csrpower_t b = std::pow(i*delta*h_r,-4./3.);
            for (uint32_t p=1; p<=maxp;p+=2) {
                csrpower_t u = std::pow(M_PI*p,2)/std::pow(2,2./3.)*b;
                Z += boost::math::airy_ai_prime(u)
                        *(boost::math::airy_ai_prime(u)
                          -j*boost::math::airy_bi_prime(u))
                  +u*boost::math::airy_ai(u)
                        *(boost::math::airy_ai(u)
                          -j*boost::math::airy_bi(u));
            }
            _data.push_back(csrpower_t(16)*Z*b*
                            csrpower_t(std::pow(M_PI,3)*std::pow(2,1./3.)/3e8));
        }
    }
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        data_buf = cl::Buffer(OCLH::context,
                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             _nmax*sizeof(impedance_t),_data.data());
    }
    #endif
}

vfps::Impedance::Impedance(const std::vector<vfps::impedance_t> &z,
                           double f_max) :
    _nmax(z.size()),
    _axis(Ruler<frequency_t>(_nmax,0,f_max,1)),
    _data(z)
{
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        data_buf = cl::Buffer(OCLH::context,
                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             _nmax*sizeof(impedance_t),_data.data());
    }
    #endif
}

vfps::Impedance::Impedance(std::string datafile, double f_max) :
    Impedance(readData(datafile),f_max)
{
}

std::vector<vfps::impedance_t> vfps::Impedance::readData(std::string fname)
{
    std::vector<vfps::impedance_t> rv;
    std::ifstream is(fname);
    size_t lineno;
    double real;
    double imag;

    while(is.good()) {
        is >> lineno >> real >> imag;
        rv.push_back(impedance_t(real,imag));
    }
    return rv;
}

uint64_t vfps::Impedance::upper_power_of_two(uint64_t v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    return v;
}
