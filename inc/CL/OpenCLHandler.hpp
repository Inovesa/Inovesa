#ifndef OPENCLHANDLER_HPP
#define OPENCLHANDLER_HPP

#define __CL_ENABLE_EXCEPTIONS

#ifdef __APPLE__
#include "inc/CL/cl.hpp"
#else
#include <CL/cl.hpp>
#endif

#include <iostream>

bool prepareCLEnvironment();

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
};

#endif // OPENCLHANDLER_HPP
