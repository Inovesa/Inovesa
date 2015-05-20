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

#include <string>

#include "fixed_point.h"

#define INOVESA_VERSION_RELEASE	0
#define INOVESA_VERSION_MINOR	6
#define INOVESA_VERSION_FIX		1

//#define INOVESA_SYNC_CL

/**
 * possible choices are:
 * 1: sum
 * 2: simpson
 */
#define INTEGRAL_TYPE 2

/**
  * possible choices are:
  * 1: no interpolation
  * 2: linear interpolation
  * 3: quadratic interpolation
  * 4: cubic interpolation
  */
#define INTERPOL_TYPE 4

/**
  * possible choices are:
  * 0: no saturation
  * 1: crop at maximum neigbouring value
  */
#define INTERPOL_SATURATING 1

/**
 * possible choices are:
 * 1: rotate on mesh
 * 2: normalized space between 0 and 1
 * 3: normalized space between -1 and +1
 */
#define ROTATION_TYPE 2

/**
  * possible choices are:
  * 1:	single-sided (only for first derivative),
  *		will "fall back" to 2 for second derivative
  * 2: two-sided (based on quadratic interpolation)
  * 3: based on cubic interpolation
  */
#define DERIVATION_TYPE 3

namespace vfps
{
#define FXP_FRACPART 54
#if FXP_FRACPART < 31
typedef fpml::fixed_point<int32_t,31-FXP_FRACPART,FXP_FRACPART> fixp32;
#endif
typedef fpml::fixed_point<int64_t,63-FXP_FRACPART,FXP_FRACPART> fixp64;

// check OpenCL code when changing (OpenCL set to uint32_t)
typedef uint32_t meshindex_t;

typedef float meshaxis_t;
typedef float meshdata_t;
typedef float interpol_t;
typedef float integral_t;
}

#endif // DEFINES_HPP
