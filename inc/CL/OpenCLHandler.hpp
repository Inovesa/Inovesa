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

#ifndef OPENCLHANDLER_HPP
#define OPENCLHANDLER_HPP
#ifdef INOVESA_USE_CL

enum class clCopyDirection {
    cpu2dev,
    dev2cpu
};

#include <GL/glew.h>

#define __CL_ENABLE_EXCEPTIONS

#if defined(__APPLE__) || defined(__MACOSX)
#include "CL/local_cl.hpp"
#else
#include <CL/cl.hpp>
#endif

#include <climits>
#include <iostream>

#include "IO/Display.hpp"


/**
 * @brief prepareCLEnvironment
 * @return true on successful initialization
 *
 * Picks the last available platform.
 * If several platforms are available,
 * AMD plattforms are prefered.
 *
 * @todo: check if OpenCL 1.0 or 2.x do the job
 * (currently "OpenCL 1." is set to be the required version)
 */
class OCLH
{
public:
	static void prepareCLEnvironment();

	static void prepareCLDevice(unsigned int device);

	static cl::Program prepareCLProg(std::string);

	static void listCLDevices();

	static bool active;

	static VECTOR_CLASS<cl::Platform> platforms;

	static cl::Context context;

	static VECTOR_CLASS<cl::Device> devices;

	/**
	 * @brief command queue for OpenCL
	 */
	static cl::CommandQueue queue;

	static bool ogl_sharing;
};

#endif // INOVESA_USE_CL
#endif // OPENCLHANDLER_HPP
