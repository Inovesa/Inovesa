#include "CL/OpenCLHandler.hpp"

bool prepareCLEnvironment()
{
	if (cl::Platform::get(&OCLH::platforms) != CL_SUCCESS) {
		return false;
	}

	cl_context_properties properties[] =
		{ CL_CONTEXT_PLATFORM, (cl_context_properties)(OCLH::platforms[0])(), 0};

	OCLH::context = cl::Context(CL_DEVICE_TYPE_ALL, properties);

	OCLH::devices = OCLH::context.getInfo<CL_CONTEXT_DEVICES>();

	OCLH::queue = cl::CommandQueue(OCLH::context, OCLH::devices[0]);

	return true;
}

VECTOR_CLASS<cl::Platform> OCLH::platforms;
cl::Context OCLH::context;
VECTOR_CLASS<cl::Device> OCLH::devices;
cl::CommandQueue OCLH::queue;
