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

#define INOVESA_RELEASE			0
#define INOVESA_VERSION_MINOR	5
#define INOVESA_VERSION_FIX		0

#include "fixed_point.h"

#define INOVESA_USE_GUI
#define INOVESA_USE_CL
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
 * 0: load "start.png"
 * 1: a square
 * 2: a 2D gaussian
 * 3: a rectengle (half)
 * 4: quarters with different patterns
 */
#define TEST_PATTERN 0

namespace vfps
{

constexpr double rotations = 1;
constexpr unsigned int patterndim_x = 512;
constexpr unsigned int patterndim_y = 4048;
constexpr unsigned int pattern_margin = 128;

constexpr unsigned int steps = 4000;

constexpr double f_s = 8.5e3;
constexpr double t_d = 0.01;

/**
 * @brief ps_xsize horizontal size of the phase space (in mesh points)
 */
constexpr unsigned int ps_xsize = 512;

/**
 * @brief ps_ysize vertical size of the phase space (in mesh points)
 */
constexpr unsigned int ps_ysize = 512;

typedef fpml::fixed_point<int32_t,2,29> fixp32;
typedef fpml::fixed_point<int64_t,34,29> fixp64;

typedef float meshaxis_t;
typedef fixp32 meshdata_t;
typedef fixp32 interpol_t;
typedef fixp64 integral_t;
}

#endif // DEFINES_HPP
