/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2014-2015: Patrik Schönfeldt                                 */
/*                                                                            */
/* This file is part of Inovesa.                                              */
/* Inovesa is free software: you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* Inovesa is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           */
/******************************************************************************/

#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <complex>
#include <string>

#include "fixed_point.h"

#define INOVESA_VERSION_RELEASE	0
#define INOVESA_VERSION_MINOR	6
#define INOVESA_VERSION_FIX		2

//#define INOVESA_SYNC_CL

/**
 * possible choices are:
 * 1: sum
 * 2: simpson
 */
#define INTEGRAL_TYPE 2

#define ROTMAP_SIZE 2

namespace vfps
{
#define FXP_FRACPART 28
#if FXP_FRACPART < 31
typedef fpml::fixed_point<int32_t,31-FXP_FRACPART,FXP_FRACPART> fixp32;
#endif
typedef fpml::fixed_point<int64_t,63-FXP_FRACPART,FXP_FRACPART> fixp64;

// has to be uint32_t (same as cl_uint) for OpenCL support
typedef uint32_t meshindex_t;

typedef double csrpower_t;
typedef std::complex<csrpower_t> impedance_t;

typedef double meshaxis_t;
typedef double meshdata_t;
typedef double interpol_t;
typedef double integral_t;
}

#endif // DEFINES_HPP