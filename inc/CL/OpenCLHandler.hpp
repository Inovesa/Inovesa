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
