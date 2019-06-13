// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include "Z/Impedance.hpp"

#include <fstream>
#include <iostream>

vfps::Impedance::Impedance(const Impedance &other)
  : Impedance( Ruler<frequency_t>(other._axis),other._data, other._oclh)
{
}

vfps::Impedance::Impedance(Ruler<frequency_t> &&axis
                          , const std::vector<vfps::impedance_t> &z
                          , oclhptr_t oclh
                          )
  : _nfreqs(z.size())
  , _axis(std::move(axis))
  , _data(z)
  , _oclh(oclh)
{
    #if INOVESA_USE_OPENCL == 1
    syncCLMem();
    #endif // INOVESA_USE_OPENCL
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

vfps::Impedance &vfps::Impedance::operator=(vfps::Impedance other)
{
    other.swap(*this);
    return *this;
}

vfps::Impedance &vfps::Impedance::operator+=(const vfps::Impedance &rhs)
{
    for (size_t i=0; i<_nfreqs; i++) {
        _data[i] += rhs._data[i];
    }
    #if INOVESA_USE_OPENCL == 1
    syncCLMem();
    #endif
    return *this;
}

void vfps::Impedance::swap(vfps::Impedance &other)
{
    std::swap(_data, other._data);
}

#if INOVESA_USE_OPENCL == 1
void vfps::Impedance::syncCLMem()
{
    if (_oclh) {
        data_buf = cl::Buffer(_oclh->context,
                              CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             _nfreqs*sizeof(impedance_t),_data.data());
    }
}
#endif

std::vector<vfps::impedance_t> vfps::Impedance::readData(std::string fname)
{
    std::vector<vfps::impedance_t> rv;
    std::ifstream is(fname);
    size_t lineno;
    size_t old_lineno(std::numeric_limits<size_t>::max());
    frequency_t real;
    frequency_t imag;

    while(is.good()) {
        is >> lineno >> real >> imag;
        if (lineno != old_lineno) {
            rv.push_back(impedance_t(real,imag));
        }
        old_lineno = lineno;
    }
    return rv;
}

constexpr double vfps::Impedance::factor4Ohms;
