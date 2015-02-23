/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2007-2009: Peter Schregle (Fixed Point Math Library)         */
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

#ifndef OPENCLHANDLER_HPP
#define OPENCLHANDLER_HPP

#define __CL_ENABLE_EXCEPTIONS

#include "CL/cl.hpp"

#include <climits>
#include <iostream>

/**
 * @brief prepareCLEnvironment
 * @param device
 * @return
 *
 * Picks the last available platform.
 * If several platforms are available,
 * AMD plattforms are prefered.
 *
 * @todo: check if OpenCL 1.0 or 2.x do the job
 * (currently "OpenCL 1." is set to be the required version)
 */
bool prepareCLEnvironment(unsigned int device=0);

class OCLH
{
public:
	static VECTOR_CLASS<cl::Platform> platforms;

	static cl::Context context;

	static VECTOR_CLASS<cl::Device> devices;

	/**
	 * @brief command queue for OpenCL
	 */
	static cl::CommandQueue queue;

	static bool ogl_sharing;
};

#endif // OPENCLHANDLER_HPP
