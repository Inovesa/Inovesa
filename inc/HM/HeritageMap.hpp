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
	HeritageMap(PhaseSpace* in, PhaseSpace* out,
				unsigned int xsize, unsigned int ysize);

	~HeritageMap();

	/**
	 * @brief apply
	 *
	 * @todo get rid of copying from/to host RAM every step
	 */
	virtual void apply();

protected:
	std::array<hi,INTERPOL_TYPE*INTERPOL_TYPE>** _heritage_map;
	std::array<hi,INTERPOL_TYPE*INTERPOL_TYPE>* const _heritage_map1D;

	const unsigned int _size;
	const unsigned int _xsize;
	const unsigned int _ysize;

	#ifdef INOVESA_USE_CL
	cl::Buffer _heritage_map1D_buf;
	cl::Kernel applyHM;

	/**
	 * @brief __initOpenCL initialize OpenCL
	 * (use after _heritage_map1D has been filled)
	 */
	void __initOpenCL();
	#endif // INOVESA_USE_CL

	PhaseSpace* _in;
	PhaseSpace* _out;
};

}

#endif // HERITAGEMAP_HPP
