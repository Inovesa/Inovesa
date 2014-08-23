#ifndef OPENCLHANDLER_HPP
#define OPENCLHANDLER_HPP

#define __CL_ENABLE_EXCEPTIONS

#include "CL/cl.hpp"

#include <iostream>

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
};

#endif // OPENCLHANDLER_HPP
