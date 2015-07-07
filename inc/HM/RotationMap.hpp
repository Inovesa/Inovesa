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

#ifndef ROTATIONMAP_HPP
#define ROTATIONMAP_HPP

#include <array>

#include "defines.hpp"
#include "HeritageMap.hpp"

namespace vfps
{

class RotationMap : public HeritageMap
{
public:
	enum class RotationCoordinates : uint_fast8_t {
		mesh = 1, // rotate on mesh
		norm_0_1 = 2, // normalized space between 0 and 1
		norm_pm1 = 3 // normalized space between -1 and +1
	};

public:
	RotationMap(PhaseSpace* in, PhaseSpace* out,
				const meshindex_t xsize, const meshindex_t ysize,
				const meshaxis_t angle, const InterpolationType it,
				const RotationCoordinates rt, bool interpol_saturating);

	/**
	 * @brief apply
	 *
	 * Saturation only makes sense with quadratic/cubic interpolation.
	 * Enabling it with linear (or without) currently crashes the program.
	 */
	void apply();

private:
	const bool _sat;

	void genHInfo(meshindex_t q_i, meshindex_t p_i, hi* myhinfo);

	const RotationCoordinates _rt;
	const meshaxis_t _cos_dt;
	const meshaxis_t _sin_dt;

	#ifdef INOVESA_USE_CL
	void genCode4HM4_2sat();

	void genCode4HM4sat();

	void genCode4Rotation();

	cl_int2 imgsize;

	cl_double2 rot;
	#endif
};

}

#endif // ROTATIONMAP_HPP
