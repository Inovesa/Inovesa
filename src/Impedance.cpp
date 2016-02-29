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
#include <fstream>

vfps::Impedance::Impedance(vfps::Impedance::ImpedanceModel model, size_t n,
                           const frequency_t f_rev, const frequency_t f_max,
                           const double h) :
    _nfreqs(n),
    _axis(Ruler<frequency_t>(_nfreqs,0,f_max,1))
{
    _data.reserve(_nfreqs);

    // according to Eq. 6.18 in Part. Acc. Vol 57, p 35 (Murpy et al.)
    constexpr impedance_t Z0 = impedance_t(306.3,176.9);

    // frequency resolution: impedance will be sampled at multiples of delta
    const frequency_t delta = f_max/f_rev/(_nfreqs-1.0);

    switch (model) {
    case ImpedanceModel::ConstantOne:
        for (size_t i=0; i<_nfreqs; i++) {
            _data.push_back(impedance_t(1,0));
        }
        break;
    case ImpedanceModel::FreeSpaceCSR:
        for (size_t i=0; i<_nfreqs; i++) {
            _data.push_back(Z0*std::pow(i*delta,csrpower_t(1.0/3.0)));
        }
        break;
    case ImpedanceModel::ParallelPlates:
        const double r_bend = physcons::c/(2*M_PI*f_rev);
        constexpr std::complex<double> j(0,1);
        for (size_t i=0; i<_nfreqs; i++) {
            std::complex<double> Z=0;
            double n = i*delta;
            double m = n*std::pow(h/r_bend,3./2.);
            const uint32_t maxp = 2*m*std::pow(r_bend/h,3./2.)*f_rev*h/physcons::c;
            double b = std::pow(m,-4./3.);
            std::complex<double> zinc=1;
            for (uint32_t p=1; p<=maxp;p+=2) {
                double u = std::pow(M_PI*p,2)/std::pow(2,2./3.)*b;
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
            Z *= 4.0*std::pow(M_PI,2)*std::pow(2,1./3.)/(physcons::epsilon0*physcons::c);
            _data.push_back(impedance_t(Z));
        }
    }
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        data_buf = cl::Buffer(OCLH::context,
                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             _nfreqs*sizeof(impedance_t),_data.data());
    }
    #endif
}

vfps::Impedance::Impedance(const std::vector<vfps::impedance_t> &z,
                           double f_max) :
    _nfreqs(z.size()),
    _axis(Ruler<frequency_t>(_nfreqs,0,f_max,1)),
    _data(z)
{
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        data_buf = cl::Buffer(OCLH::context,
                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             _nfreqs*sizeof(impedance_t),_data.data());
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
