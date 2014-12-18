#include "CL/OpenCLHandler.hpp"

bool prepareCLEnvironment(unsigned int device)
{
	if (cl::Platform::get(&OCLH::platforms) != CL_SUCCESS) {
		return false;
	}

    std::string available_clversion;
    std::string needed_clversion = "OpenCL 1.2";

    unsigned int plati = 0;
    std::string cl_plat_vendor;
    for (unsigned int p=0; p<OCLH::platforms.size(); p++) {
        OCLH::platforms[p].getInfo(CL_PLATFORM_VERSION,&available_clversion);
        if (available_clversion.substr(0,needed_clversion.length()) == needed_clversion) {
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

	std::string devicename;
	if (device == UINT_MAX) {
		std::cout << "Available OpenCL devices:" << std::endl;
        for (unsigned int d=0; d<OCLH::devices.size(); d++) {
            OCLH::devices[d].getInfo(CL_DEVICE_NAME,&devicename);
            std::cout << "==Device " << d+1 << "==" << std::endl
					  << devicename << " ("
                      << OCLH::devices[d].getInfo<CL_DEVICE_EXTENSIONS>()
					  << ")" << std::endl;
		}
		std::cout << "Choose which device to use: ";
		std::cin >> device;
		device = (device-1)%OCLH::devices.size();
	}
    OCLH::devices[device].getInfo<std::string>(CL_DEVICE_NAME,&devicename);
	OCLH::queue = cl::CommandQueue(OCLH::context, OCLH::devices[device]);
	// cl_VENDOR_gl_sharing is present, when string contains the substring
	OCLH::ogl_sharing
			= OCLH::devices[device].getInfo<CL_DEVICE_EXTENSIONS>().find(
				"_gl_sharing") != std::string::npos;
    #ifdef DEBUG
    std::cout << "Available OpenCL Devices:" << std::endl;
    for (unsigned int p=0; p<OCLH::platforms.size(); p++) {
        cl_context_properties tmp_properties[] =
            { CL_CONTEXT_PLATFORM, (cl_context_properties)(OCLH::platforms[p])(), 0};
        cl::Context tmp_context = cl::Context(CL_DEVICE_TYPE_ALL, tmp_properties);
        VECTOR_CLASS<cl::Device> tmp_devices = tmp_context.getInfo<CL_CONTEXT_DEVICES>();
        for (unsigned int d=0; d<OCLH::devices.size(); d++) {
            if (p == plati && d == device) {
                std::cout << "(*) ";
            } else {
                std::cout << "( ) ";
            }
            std::cout << tmp_devices[d].getInfo<CL_DEVICE_NAME>()
                      << " (on \"" <<  OCLH::platforms[p].getInfo<CL_PLATFORM_NAME>() << "\")"
                      << " by \"" <<  OCLH::platforms[p].getInfo<CL_PLATFORM_VENDOR>() << '\"'
                      << " offering \"" << tmp_devices[d].getInfo<CL_DEVICE_VERSION>() << '\"'
//                      << " with " << tmp_devices[d].getInfo<CL_DEVICE_EXTENSIONS>()
                      << std::endl;
        }
    }
    #else
    std::cout << "Using " << devicename << " for OpenCL." << std::endl;
    #endif

	return true;
}

VECTOR_CLASS<cl::Platform> OCLH::platforms;
cl::Context OCLH::context;
VECTOR_CLASS<cl::Device> OCLH::devices;
cl::CommandQueue OCLH::queue;
bool OCLH::ogl_sharing;
