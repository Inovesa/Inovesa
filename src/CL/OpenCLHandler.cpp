#include "CL/OpenCLHandler.hpp"

bool prepareCLEnvironment(unsigned int device)
{
	if (cl::Platform::get(&OCLH::platforms) != CL_SUCCESS) {
		return false;
	}

	cl_context_properties properties[] =
		{ CL_CONTEXT_PLATFORM, (cl_context_properties)(OCLH::platforms[0])(), 0};

	OCLH::context = cl::Context(CL_DEVICE_TYPE_ALL, properties);

	OCLH::devices = OCLH::context.getInfo<CL_CONTEXT_DEVICES>();

	std::string devicename;
	if (device == UINT_MAX) {
		std::cout << "Available OpenCL devices:" << std::endl;
		for (unsigned int i=0; i<OCLH::devices.size(); i++) {
			OCLH::devices[i].getInfo<std::string>(CL_DEVICE_NAME,&devicename);
			std::cout << "==Device " << i+1 << "==" << std::endl
					  << devicename << " ("
					  << OCLH::devices[i].getInfo<CL_DEVICE_EXTENSIONS>()
					  << ")" << std::endl;
		}
		std::cout << "Choose which device to use: ";
		std::cin >> device;
		device = (device-1)%OCLH::devices.size();
	}
	OCLH::devices[device].getInfo<std::string>(CL_DEVICE_NAME,&devicename);
	std::cout << "Using " << devicename << " for OpenCL." << std::endl;
	OCLH::queue = cl::CommandQueue(OCLH::context, OCLH::devices[device]);
	// cl_VENDOR_gl_sharing is present, when string contains the substring
	OCLH::ogl_sharing
			= OCLH::devices[device].getInfo<CL_DEVICE_EXTENSIONS>().find(
				"_gl_sharing") != std::string::npos;

	return true;
}

VECTOR_CLASS<cl::Platform> OCLH::platforms;
cl::Context OCLH::context;
VECTOR_CLASS<cl::Device> OCLH::devices;
cl::CommandQueue OCLH::queue;
bool OCLH::ogl_sharing;
