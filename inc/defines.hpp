/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2017: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2017: Karlsruhe Institute of Technology                 *
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

#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <complex>
#include <string>
#include <fftw3.h>

#include "fixed_point.h"

#define INOVESA_VERSION_RELEASE 0
#define INOVESA_VERSION_MINOR   14

/* pre release:
 * -1: alpha
 * -2: beta
 * -3 RC1
 * -n RC(n-2)
 */
#define INOVESA_VERSION_FIX     -11

//#define INOVESA_SYNC_CL

namespace vfps
{
#define FXP_FRACPART 28
#if FXP_FRACPART < 31
typedef fpml::fixed_point<int32_t,31-FXP_FRACPART,FXP_FRACPART> fixp32;
#endif
typedef fpml::fixed_point<int64_t,63-FXP_FRACPART,FXP_FRACPART> fixp64;

// has to be uint32_t (same as cl_uint) for OpenCL support
typedef uint32_t meshindex_t;

/* currently all of the below types have to be the same
 * use this to switch them all */
typedef float data_t;

typedef data_t csrpower_t;
typedef std::complex<csrpower_t> impedance_t;

typedef data_t frequency_t;
typedef data_t meshaxis_t;
typedef data_t meshdata_t;
typedef data_t interpol_t;
typedef data_t integral_t;
typedef data_t timeaxis_t;

typedef integral_t projection_t;

inline bool isOfFileType(std::string ending, std::string fname)
{
    return ( fname.size() > ending.size() &&
        std::equal(ending.rbegin(), ending.rend(),fname.rbegin()));
}

namespace physcons
{
/// speed of light (in m/s)
constexpr double c=2.99792458e8;

/// Vacuum permittivity (in F/m)
constexpr double epsilon0 = 8.854187817e-12;

/// charge of one electron (in C)
constexpr double e=1.602e-19;

/// Alfven current (in A)
constexpr double IAlfven = 17045;

/// electron rest energy (in eV)
constexpr double me = 510998.9;

/// magnetic constant (in H/m)
constexpr double mu0 = 1/(epsilon0*c*c);
} // namespace physcons

} // namespace vfps

#endif // DEFINES_HPP
