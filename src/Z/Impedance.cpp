/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2016: Karlsruhe Institute of Technology                 *
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

vfps::Impedance::Impedance(const std::vector<vfps::impedance_t> &z,
                           const frequency_t f_max) :
    _nfreqs(z.size()),
    _axis(Ruler<frequency_t>(_nfreqs,0,f_max,1)),
    _data(z)
{
    syncCLMem();
}

vfps::Impedance::Impedance(std::string datafile, double f_max) :
    Impedance(readData(datafile),f_max)
{
}

void vfps::Impedance::syncCLMem()
{
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        data_buf = cl::Buffer(OCLH::context,
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
