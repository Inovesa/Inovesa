/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2015: Patrik Schönfeldt                                 *
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

vfps::Impedance::Impedance(vfps::Impedance::ImpedanceModel model, size_t n,
                           double f_rev, double f_max, bool roundup) :
    _nmax(roundup? upper_power_of_two(n) : n),
    _axis(Ruler<frequency_t>(_nmax,0,f_max,1))
{
    _data.reserve(_nmax);

    // according to Eq. 6.18 in Part. Acc. Vol 57, p 35 (Murpy et al.)
    constexpr impedance_t freespacecoeff = impedance_t(250.1,176.9);

    // frequency resolution: impedance will be sampled at multiples of deltaf
    const frequency_t deltaf = f_max/f_rev*(_nmax/(_nmax-1.0));

    switch (model) {
    case ImpedanceModel::FreeSpace:
        for (size_t i=0; i<_nmax; i++) {
            _data.push_back(freespacecoeff*std::pow(i*deltaf,
                                                    csrpower_t(1.0/3.0)));
        }
    break;
    }
}

vfps::Impedance::Impedance(const std::vector<vfps::impedance_t> &z,
                           double f_max) :
    _nmax(z.size()),
    _axis(Ruler<frequency_t>(_nmax,0,f_max,1)),
    _data(z)
{
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
