// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#include "Z/Impedance.hpp"

#include <fstream>

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
