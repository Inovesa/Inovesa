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

#include "PhaseSpace.hpp"

namespace vfps
{

class HeritageMap
{
protected:
	typedef struct {
		unsigned int index;
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
				unsigned int xsize, unsigned int ysize,
				unsigned int interpoints);

	~HeritageMap();

	/**
	 * @brief apply
	 */
	virtual void apply();

protected:
	/**
	 * @brief _ip holds the number of points used for interpolation
	 */
	const unsigned int _ip;

	hi*** _heritage_map;
	hi** const _heritage_map1D;
	hi* const _hinfo;

	const unsigned int _size;
	const unsigned int _xsize;
	const unsigned int _ysize;

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
