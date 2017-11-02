/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2012-2016: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2016: Karlsruhe Institute of Technology                 *
 *                                                                            *
 * This file is part of Inovesa.                                              *
 * Inovesa is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Inovesa is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           *
 ******************************************************************************/

#ifndef OPENCLHANDLER_HPP
#define OPENCLHANDLER_HPP
#ifdef INOVESA_USE_CL

enum class clCopyDirection {
    cpu2dev,
    dev2cpu
};

#include <GL/glew.h>

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 110
#define CL_HPP_TARGET_OPENCL_VERSION 120

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include "CL/local_cl.hpp"
#pragma GCC diagnostic pop

#define INOVESA_ENABLE_CLPROFILING
#ifdef INOVESA_ENABLE_CLPROFILING
#include "CL/CLProfiler.hpp"
#endif

#ifdef INOVESA_USE_CLFFT
#include <clFFT.h>
#endif // INOVESA_USE_CLFFT

#include <climits>
#include <list>
#include <iostream>

#define INOVESA_ENABLE_CLPROFILING

#ifdef INOVESA_ENABLE_CLPROFILING
#include "CL/CLProfiler.hpp"
#endif // INOVESA_ENABLE_CLPROFILING

/**
 * Picks the last available platform.
 * If several platforms are available,
 * AMD plattforms are prefered.
 *
 * @todo: update to recent OpenCL version
 *        (currently OpenCL 1.1 or 1.2 is  required)
 */
class OCLH
{
public:
    static void prepareCLEnvironment(bool glsharing, uint32_t device);

    static cl::Program prepareCLProg(std::string);

    static void teardownCLEnvironment();

    static void listCLDevices();

    static bool active;

    static cl::vector<cl::Platform> platforms;

    static cl::Context context;

    static cl::vector<cl::Device> devices;

    static cl_device_type devicetype;

    /**
     * @brief command queue for OpenCL
     */
    static cl::CommandQueue queue;

    static bool ogl_sharing;


    #ifdef INOVESA_ENABLE_CLPROFILING
    static std::list<vfps::CLTiming> timings;

    static cl::Event init;
    #endif // INOVESA_ENABLE_CLPROFILING

public:
    /**
     * This wrapper function allows to centrally controll queuing kernels.
     * At the moment, it just forwards the arguments.
     */
    static inline void
    enqueueNDRangeKernel(const cl::Kernel& kernel,
                         const cl::NDRange& offset,
                         const cl::NDRange& global,
                         const cl::NDRange& local = cl::NullRange,
                         const cl::vector<cl::Event>* events = nullptr,
                         cl::Event* event = nullptr)
    {
        queue.enqueueNDRangeKernel(kernel,offset,global,local,events,event);
    }

    /**
     * This wrapper function allows to centrally controll queuing copyBuffer
     */
    static inline void
    enqueueCopyBuffer(const cl::Buffer& src,
                      const cl::Buffer& dst,
                      cl::size_type src_offset,
                      cl::size_type dst_offset,
                      cl::size_type size,
                      const cl::vector<cl::Event>* events = nullptr,
                      cl::Event* event = nullptr)
    {
        queue.enqueueCopyBuffer(src, dst, src_offset,dst_offset,size,
                                events, event);
    }

private:
#ifdef INOVESA_USE_CLFFT
        static clfftSetupData fft_setup;
#endif // INOVESA_USE_CLFFT

    static const std::string custom_datatypes;

    static std::string datatype_aliases();
};

#endif // INOVESA_USE_CL
#endif // OPENCLHANDLER_HPP
