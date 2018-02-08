/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2012-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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
#ifdef INOVESA_USE_OPENCL

enum class clCopyDirection {
    cpu2dev,
    dev2cpu
};

#include <GL/glew.h>

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 110
#define CL_HPP_TARGET_OPENCL_VERSION 120

#include "CL/local_cl.hpp"

#ifdef INOVESA_ENABLE_CLPROFILING
#include "CL/CLProfiler.hpp"
#endif

#ifdef INOVESA_USE_CLFFT
#include <clFFT.h>
#endif // INOVESA_USE_CLFFT

#include <climits>
#include <list>
#include <iostream>

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
    /**
     * @brief prepareCLEnvironment
     * @param device
     * @param glsharing needs to be implemented
     */
    static void prepareCLEnvironment( uint32_t device
                                    #ifdef INOVESA_USE_OPENGL
                                    , bool glsharing
                                    #endif
                                    );

    static cl::Program prepareCLProg(std::string);

    static void teardownCLEnvironment();

    static void teardownCLEnvironment(cl::Error& e);

    static void listCLDevices();

    static bool active;

    static cl::Context context;

private:
    static cl::vector<cl::Platform> platforms;

    static cl::vector<cl::Device> devices;

    static cl_device_type devicetype;

    /**
     * @brief command queue for OpenCL
     */
    static cl::CommandQueue queue;

    #ifdef INOVESA_USE_OPENGL
    static bool ogl_sharing;
    #endif // INOVESA_USE_OPENGL

    #ifdef INOVESA_ENABLE_CLPROFILING
    static void saveProfilingInfo(std::string fname);

    static std::list<vfps::CLTiming> timingInfo;

    static cl::Event init;
    #endif // INOVESA_ENABLE_CLPROFILING

public:
    #ifdef INOVESA_USE_CLFFT
    static inline void
    bakeClfftPlan(clfftPlanHandle plHandle)
    {
        clfftBakePlan(plHandle,1,&queue(), nullptr, nullptr);
    }

    static inline void
    enqueueDFT(clfftPlanHandle plHandle,
               clfftDirection dir,
               cl::Buffer inputBuffer,
               cl::Buffer outputBuffer)
    {
        #ifdef INOVESA_ENABLE_CLPROFILING
        cl::Event* event = new cl::Event();
        timingsDFT.push_back(event);
        #endif // INOVESA_ENABLE_CLPROFILING
        clfftEnqueueTransform(plHandle,dir,1,&queue(),0,nullptr,
                              #ifdef INOVESA_ENABLE_CLPROFILING
                              &(*event)(),
                              #else
                              nullptr,
                              #endif
                              &inputBuffer(),&outputBuffer(),nullptr);
    }
    #endif // INOVESA_USE_CLFFT

    static inline void
    enqueueBarrier()
    {
        #ifdef CL_VERSION_1_2
        OCLH::queue.enqueueBarrierWithWaitList();
        #else // CL_VERSION_1_2
        OCLH::queue.enqueueBarrier();
        #endif // CL_VERSION_1_2
    }

    /**
     * This wrapper function allows to centrally controll queuing kernels.
     */
    static inline void
    enqueueNDRangeKernel(const cl::Kernel& kernel,
                         const cl::NDRange& offset,
                         const cl::NDRange& global,
                         const cl::NDRange& local = cl::NullRange,
                         const cl::vector<cl::Event>* events = nullptr,
                         cl::Event* event = nullptr
                         #ifdef INOVESA_ENABLE_CLPROFILING
                         , cl::vector<cl::Event*>* timings = nullptr
                         #endif // INOVESA_ENABLE_CLPROFILING
                         )
    {
        #ifdef INOVESA_ENABLE_CLPROFILING
        if (event == nullptr) {
            event = new cl::Event();
        }
        if (timings == nullptr) {
            timingsExecute.push_back(event);
        } else {
            timings->push_back(event);
        } // INOVESA_ENABLE_CLPROFILING
        #endif
        queue.enqueueNDRangeKernel(kernel,offset,global,local,events,event);

        enqueueBarrier();
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
                      cl::Event* event = nullptr
                      #ifdef INOVESA_ENABLE_CLPROFILING
                      , cl::vector<cl::Event*>* timings = nullptr
                      #endif // INOVESA_ENABLE_CLPROFILING
                      )
    {
        #ifdef INOVESA_ENABLE_CLPROFILING
        if (event == nullptr) {
            event = new cl::Event();
        }
        if (timings == nullptr) {
            timingsCopy.push_back(event);
        } else {
            timings->push_back(event);
        }
        #endif // INOVESA_ENABLE_CLPROFILING
        queue.enqueueCopyBuffer(src, dst, src_offset,dst_offset,size,
                                events, event);
    }

    /**
     * This wrapper function allows to centrally controll queuing readBuffer
     */
    static inline void
    enqueueReadBuffer(const cl::Buffer& buffer,
                       cl_bool blocking,
                       cl::size_type src_offset,
                       cl::size_type size,
                       void* ptr,
                       const cl::vector<cl::Event>* events = nullptr,
                       cl::Event* event = nullptr
                     #ifdef INOVESA_ENABLE_CLPROFILING
                     , cl::vector<cl::Event*>* timings = nullptr
                     #endif // INOVESA_ENABLE_CLPROFILING
                     )
    {
        #ifdef INOVESA_ENABLE_CLPROFILING
        if (event == nullptr) {
            event = new cl::Event();
        }
        if (timings == nullptr) {
            timingsRead.push_back(event);
        } else {
            timings->push_back(event);
        }
        #endif // INOVESA_ENABLE_CLPROFILING
        queue.enqueueReadBuffer(buffer, blocking, src_offset,size,ptr,
                                events, event);
    }

    /**
     * This wrapper function allows to centrally controll queuing writeBuffer
     */
    static inline void
    enqueueWriteBuffer( const cl::Buffer& buffer
                      , cl_bool blocking
                      , cl::size_type src_offset
                      , cl::size_type size
                      , const void* ptr
                      , const cl::vector<cl::Event>* events = nullptr
                      , cl::Event* event = nullptr
                      #ifdef INOVESA_ENABLE_CLPROFILING
                      , cl::vector<cl::Event*>* timings = nullptr
                      #endif // INOVESA_ENABLE_CLPROFILING
                      )
    {
        #ifdef INOVESA_ENABLE_CLPROFILING
        if (event == nullptr) {
            event = new cl::Event();
        }
        if (timings == nullptr) {
            timingsWrite.push_back(event);
        } else {
            timings->push_back(event);
        }
        #endif // INOVESA_ENABLE_CLPROFILING
        queue.enqueueWriteBuffer(buffer, blocking, src_offset,size,ptr,
                                events, event);
    }

    static inline void finish()
    {
        OCLH::queue.finish();
    }

    static inline void flush()
    {
        OCLH::queue.flush();
    }

    #ifdef INOVESA_ENABLE_CLPROFILING
    static void saveTimings(cl::vector<cl::Event *> *evts, std::string name);
    #endif // INOVESA_ENABLE_CLPROFILING

private:
    #ifdef INOVESA_USE_CLFFT
    static clfftSetupData fft_setup;
    #endif // INOVESA_USE_CLFFT

    #ifdef INOVESA_ENABLE_CLPROFILING
    static cl::vector<cl::Event*> timingsCopy;
    static cl::vector<cl::Event*> timingsDFT;
    static cl::vector<cl::Event*> timingsExecute;
    static cl::vector<cl::Event*> timingsRead;
    static cl::vector<cl::Event*> timingsWrite;
    #endif // INOVESA_ENABLE_CLPROFILING

    static const std::string custom_datatypes;

    static std::string datatype_aliases();
};

#endif // INOVESA_USE_OPENCL
#endif // OPENCLHANDLER_HPP
