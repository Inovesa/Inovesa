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
#include "IO/Display.hpp"

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

    clfftInitSetupData(&fft_setup);
    clfftSetup(&fft_setup);

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
    code = datatype_aliases()+custom_datatypes+code;
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

void OCLH::teardownCLEnvironment()
{
    clfftTeardown();
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

clfftSetupData OCLH::fft_setup;

const std::string OCLH::custom_datatypes = R"(
        typedef struct { data_t real; data_t imag; } impedance_t;
        inline impedance_t cmult(const impedance_t a, const impedance_t b) {
            impedance_t rv;
            rv.real = a.real*b.real - a.imag*b.imag;
            rv.imag = a.imag*b.real + a.real*b.imag;
            return rv;
        }

        )";

std::string OCLH::datatype_aliases()
{
    std::string code;
    if (std::is_same<vfps::meshdata_t,float>::value) {
    code +=
        "typedef float data_t;\n"
        "typedef float2 data2_t;\n"
        "typedef float3 data3_t;\n"
        "typedef float4 data4_t;\n"
        "float mult(float x, float y);"
        "float mult(float x, float y) { return x*y; }\n";
    } else  if (std::is_same<vfps::meshdata_t,double>::value) {
    code +=
        "typedef double data_t;\n"
        "typedef double2 data2_t;\n"
        "typedef double3 data3_t;\n"
        "typedef double4 data4_t;\n"
        "double mult(double x, double y);"
        "double mult(double x, double y) { return x*y; }\n";
    } else {
        std::stringstream fxp_fracpart;
        fxp_fracpart << FXP_FRACPART;

        code +=    "__constant int fracpart="+fxp_fracpart.str()+";\n";
        #if FXP_FRACPART < 31
        if (std::is_same<vfps::meshdata_t,vfps::fixp32>::value) {
        code +=
            "typedef int data_t;\n"
            "typedef int2 data2_t;\n"
            "typedef int3 data3_t;\n"
            "typedef int4 data4_t;\n"
            "int mult(int x, int y){return ((long)(x)*(long)(y))>>fracpart;}\n";
        } else
        #endif
        if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value) {
        code +=
            "typedef long data_t;\n"
            "typedef long2 data2_t;\n"
            "typedef long3 data3_t;\n"
            "typedef long4 data4_t;\n"
            "long mult(long x, long y) {"
            "return ((mul_hi(x,y) << (64-fracpart)) + ((x*y) >> fracpart));}\n";
        }
    }
    return code;
}

#endif
