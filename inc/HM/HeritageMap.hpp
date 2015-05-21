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

#ifndef HERITAGEMAP_HPP
#define HERITAGEMAP_HPP

#include "defines.hpp"
#include "PhaseSpace.hpp"

namespace vfps
{

class HeritageMap
{
public:
	enum InterpolationType : uint_fast8_t {
		none = 1,
		linear = 2,
		quadratic = 3,
		cubic = 4
	};

protected:
	typedef struct {
		meshindex_t index;
		interpol_t weight;
	} hi;

public:
	/**
	 * @brief HeritageMap
	 * @param in
	 * @param out
	 * @param xsize
	 * @param ysize
	 * @param interpoints number of points used for interpolation
	 */
	HeritageMap(PhaseSpace* in, PhaseSpace* out,
				size_t xsize, size_t ysize,
				uint_fast8_t interpoints, uint_fast8_t intertype);

	~HeritageMap();

	/**
	 * @brief apply
	 */
	virtual void apply();

protected:
	/**
	 * @brief _ip holds the total number of points used for interpolation
	 */
	#ifdef INOVESA_USE_CL
	const cl_uint _ip;
	#else
	const uint_fast8_t _ip;
	#endif

	/**
	 * @brief _ip holds the per dimension number
	 *			of points used for interpolation
	 */
	const uint_fast8_t _it;

	/**
	 * @brief _hinfo
	 */
	hi* const _hinfo;

	/**
	 * @brief _size size of the HeritageMap (_xsize*_ysize)
	 */
	const meshindex_t _size;

	/**
	 * @brief _xsize horizontal size of the HeritageMap
	 */
	const meshindex_t _xsize;

	/**
	 * @brief _ysize vertical size of the HeritageMap
	 */
	const meshindex_t _ysize;

	#ifdef INOVESA_USE_CL
	/**
	 * @brief _hi_buf buffer for heritage information
	 */
	cl::Buffer _hi_buf;

	/**
	 * @brief applyHM
	 */
	cl::Kernel applyHM;

	#endif // INOVESA_USE_CL

	PhaseSpace* _in;
	PhaseSpace* _out;
};

}

#endif // HERITAGEMAP_HPP
