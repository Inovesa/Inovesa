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
 *along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.            *
 ******************************************************************************/

#include "Z/Impedance.hpp"

#include <fstream>

vfps::Impedance::Impedance(const vfps::Impedance &other)
  : Impedance( other._axis,other._data, other._oclh)
{
}

vfps::Impedance::Impedance( Ruler<vfps::frequency_t> axis
                          , const std::vector<vfps::impedance_t> &z
                          , oclhptr_t oclh
                          )
  : _nfreqs(z.size())
  , _axis(axis)
  , _data(z)
  , _oclh(oclh)
{
    syncCLMem();
}

vfps::Impedance::Impedance( const std::vector<vfps::impedance_t> &z
                          , const frequency_t f_max
                          , oclhptr_t oclh
                          )
  : Impedance( Ruler<frequency_t>(z.size(),0,f_max,{{"Hertz",1}}),z
             , oclh
             )
{
}

vfps::Impedance::Impedance( const size_t nfreqs
                          , const vfps::frequency_t f_max
                          , oclhptr_t oclh
                          )
  : Impedance( Ruler<frequency_t>(nfreqs,0,f_max,{{"Hertz",1}})
             , std::vector<impedance_t>(nfreqs,0)
             , oclh
             )
{
}

vfps::Impedance::Impedance( std::string datafile
                          , double f_max
                          , oclhptr_t oclh
                          )
  : Impedance( readData(datafile),f_max, oclh)
{
}

vfps::Impedance &vfps::Impedance::operator+=(const vfps::Impedance &rhs)
{
    for (size_t i=0; i<_nfreqs; i++) {
        _data[i] += rhs._data[i];
    }
    syncCLMem();
    return *this;
}

void vfps::Impedance::syncCLMem()
{
    #if INOVESA_USE_OPENCL == 1
    if (_oclh) {
        data_buf = cl::Buffer(_oclh->context,
                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             _nfreqs*sizeof(impedance_t),_data.data());
    }
    #endif
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

constexpr double vfps::Impedance::factor4Ohms;
