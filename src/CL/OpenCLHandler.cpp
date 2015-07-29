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

#ifdef INOVESA_USE_CL

#include "CL/OpenCLHandler.hpp"

void OCLH::prepareCLEnvironment()
{
	cl::Platform::get(&OCLH::platforms);

    std::string available_clversion;
	std::string needed_clversion = "OpenCL 1.";

    unsigned int plati = 0;
    std::string cl_plat_vendor;
    for (unsigned int p=0; p<OCLH::platforms.size(); p++) {
        OCLH::platforms[p].getInfo(CL_PLATFORM_VERSION,&available_clversion);
        if (available_clversion.substr(0,needed_clversion.length()) == needed_clversion) {
			// stick with AMD (if available)
			if (cl_plat_vendor != "Advanced Micro Devices, Inc.") {
                plati = p;
                cl_plat_vendor = OCLH::platforms[p].getInfo<CL_PLATFORM_VENDOR>();
			}
        }
    }

    cl_context_properties properties[] =
        { CL_CONTEXT_PLATFORM, (cl_context_properties)(OCLH::platforms[plati])(), 0};

	OCLH::context = cl::Context(CL_DEVICE_TYPE_ALL, properties);

	OCLH::devices = OCLH::context.getInfo<CL_CONTEXT_DEVICES>();
}

void
OCLH::prepareCLDevice(unsigned int device)
{
	OCLH::queue = cl::CommandQueue(OCLH::context, OCLH::devices[device]);
	// cl_VENDOR_gl_sharing is present, when string contains the substring
	OCLH::ogl_sharing
			= OCLH::devices[device].getInfo<CL_DEVICE_EXTENSIONS>().find(
				"_gl_sharing") != std::string::npos;
	std::string devicename;
	OCLH::devices[device].getInfo<std::string>(CL_DEVICE_NAME,&devicename);
	vfps::Display::printText("Initialized \""+devicename+"\" for use with OpenCL.");
}

cl::Program OCLH::prepareCLProg(std::string code)
{
	cl::Program::Sources source(1,std::make_pair(code.c_str(),code.length()));
	cl::Program p(OCLH::context, source);
	try {
		p.build(OCLH::devices);
	} catch (cl::Error &e) {
		e.what();
	#ifdef DEBUG
	// in debug builds, CL code and build log should always be displayed
	}
	#endif  // DEBUG
		std::cout	<< "===== OpenCL Code =====\n"
					<< code << std::endl;
		std::cout	<< "===== OpenCL Build Log =====\n"
					<< p.getBuildInfo<CL_PROGRAM_BUILD_LOG>(OCLH::devices[0])
					<< std::endl;
	#ifndef DEBUG
	// in release builds show CL code and build log only when error occured
	}
	#endif // DEBUG

	return p;
}

void OCLH::listCLDevices()
{
	std::cout << "Available OpenCL Devices:" << std::endl;
	for (unsigned int p=0; p<OCLH::platforms.size(); p++) {
		cl_context_properties tmp_properties[] =
			{ CL_CONTEXT_PLATFORM, (cl_context_properties)(OCLH::platforms[p])(), 0};
		cl::Context tmp_context = cl::Context(CL_DEVICE_TYPE_ALL, tmp_properties);
		VECTOR_CLASS<cl::Device> tmp_devices = tmp_context.getInfo<CL_CONTEXT_DEVICES>();
		for (unsigned int d=0; d<OCLH::devices.size(); d++) {
			std::cout	<< "=> " << d+1 << ": "
						<< tmp_devices[d].getInfo<CL_DEVICE_NAME>() << " ("
						<< tmp_devices[d].getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>()/0x100000 << " MiB, "
						<< tmp_devices[d].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()
						<< " CU on \"" <<  OCLH::platforms[p].getInfo<CL_PLATFORM_NAME>() << "\")"
						<< " by \"" <<  OCLH::platforms[p].getInfo<CL_PLATFORM_VENDOR>() << '\"' << std::endl
						<< " offering \"" << tmp_devices[d].getInfo<CL_DEVICE_VERSION>() << '\"'
						<< " with " << tmp_devices[d].getInfo<CL_DEVICE_EXTENSIONS>()
						<< std::endl;
		}
	}
}

bool OCLH::active;
VECTOR_CLASS<cl::Platform> OCLH::platforms;
cl::Context OCLH::context;
VECTOR_CLASS<cl::Device> OCLH::devices;
cl::CommandQueue OCLH::queue;
bool OCLH::ogl_sharing;

#endif
