// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include <memory>

#if INOVESA_USE_OPENGL == 1
#include <GL/glew.h>
#endif // INOVESA_USE_OPENGL

#if INOVESA_USE_OPENCL == 1
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 110
#define CL_HPP_TARGET_OPENCL_VERSION 120

#include "CL/local_cl.hpp"
#endif // INOVESA_USE_OPENCL

// shared OpenCL OpenGL pointer with automatic fallback
#if INOVESA_USE_OPENGL == 1
namespace vfps{
#if INOVESA_USE_OPENCL == 1
typedef cl_GLuint clgluint;
#else
typedef GLuint clgluint;
#endif // INOVESA_USE_OPENCL
} // namespace vfps
#endif // INOVESA_USE_OPENGL

// negation for preprocessor to have short alternative on top
#if !(INOVESA_USE_OPENCL == 1)
typedef std::nullptr_t oclhptr_t;
#else

class OCLH;
typedef std::shared_ptr<OCLH> oclhptr_t;

#if INOVESA_ENABLE_CLPROFILING == 1
#include "CL/CLProfiler.hpp"
#endif

#if INOVESA_USE_CLFFT == 1
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
    enum class clCopyDirection {
        cpu2dev,
        dev2cpu
    };

    /**
     * @brief prepareCLEnvironment
     * @param device
     * @param glsharing needs to be implemented
     */
    OCLH( uint32_t device, bool glsharing=false);

    cl::Program prepareCLProg(std::string);

    ~OCLH();

    void teardownCLEnvironment(cl::Error& e);

     static void listCLDevices();

    cl::Context context;

    /**
     * @brief OpenGLSharing
     * @return status of OpenGL sharing (always false for build w/o OpenGL)
     *
     */
    bool OpenGLSharing() const
        { return ogl_sharing; }

private:
    cl::Platform _platform;

    cl::vector<cl::Device> _devices;

    cl::Device _device;

    cl_device_type devicetype;

    /**
     * @brief command queue for OpenCL
     */
    cl::CommandQueue queue;

    bool ogl_sharing;

    #if INOVESA_ENABLE_CLPROFILING == 1
    void saveProfilingInfo(std::string fname);

    std::list<vfps::CLTiming> timingInfo;

    cl::Event init;
    #endif // INOVESA_ENABLE_CLPROFILING

public:
    #if INOVESA_USE_CLFFT == 1
    /**
     * @brief bakeClfftPlan wrapper for clfftBakePlan
     * @param plHandle
     */
    inline void
    bakeClfftPlan(clfftPlanHandle plHandle)
    {
        clfftBakePlan(plHandle,1,&queue(), nullptr, nullptr);
    }

    inline void
    enqueueDFT(clfftPlanHandle plHandle,
               clfftDirection dir,
               cl::Buffer inputBuffer,
               cl::Buffer outputBuffer)
    {
        #if INOVESA_ENABLE_CLPROFILING == 1
        cl::Event* event = new cl::Event();
        timingsDFT.push_back(event);
        #endif // INOVESA_ENABLE_CLPROFILING
        clfftEnqueueTransform(plHandle,dir,1,&queue(),0,nullptr,
                              #if INOVESA_ENABLE_CLPROFILING == 1
                              &(*event)(),
                              #else
                              nullptr,
                              #endif
                              &inputBuffer(),&outputBuffer(),nullptr);
    }
    #endif // INOVESA_USE_CLFFT

    inline void
    enqueueBarrier()
    {
        #if CL_VERSION_1_2 == 1
        queue.enqueueBarrierWithWaitList();
        #elif CL_VERSION_1_1 == 1
        queue.enqueueBarrier();
        #endif // CL_VERSION_1_2
    }

    /**
     * This wrapper function allows to centrally controll queuing kernels.
     */
    inline void
    enqueueNDRangeKernel(const cl::Kernel& kernel,
                         const cl::NDRange& offset,
                         const cl::NDRange& global,
                         const cl::NDRange& local = cl::NullRange,
                         const cl::vector<cl::Event>* events = nullptr,
                         cl::Event* event = nullptr
                         #if INOVESA_ENABLE_CLPROFILING == 1
                         , cl::vector<cl::Event*>* timings = nullptr
                         #endif // INOVESA_ENABLE_CLPROFILING
                         )
    {
        #if INOVESA_ENABLE_CLPROFILING == 1
        if (event == nullptr) {
            event = new cl::Event();
        }
        if (timings == nullptr) {
            timingsExecute.push_back(event);
        } else {
            timings->push_back(event);
        }
        #endif // INOVESA_ENABLE_CLPROFILING
        queue.enqueueNDRangeKernel(kernel,offset,global,local,events,event);

        enqueueBarrier();
    }

    /**
     * This wrapper function allows to centrally controll queuing copyBuffer
     */
    inline void
    enqueueCopyBuffer(const cl::Buffer& src,
                      const cl::Buffer& dst,
                      cl::size_type src_offset,
                      cl::size_type dst_offset,
                      cl::size_type size,
                      const cl::vector<cl::Event>* events = nullptr,
                      cl::Event* event = nullptr
                      #if INOVESA_ENABLE_CLPROFILING == 1
                      , cl::vector<cl::Event*>* timings = nullptr
                      #endif // INOVESA_ENABLE_CLPROFILING
                      )
    {
        #if INOVESA_ENABLE_CLPROFILING == 1
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
    inline void
    enqueueReadBuffer(const cl::Buffer& buffer,
                       cl_bool blocking,
                       cl::size_type src_offset,
                       cl::size_type size,
                       void* ptr,
                       const cl::vector<cl::Event>* events = nullptr,
                       cl::Event* event = nullptr
                     #if INOVESA_ENABLE_CLPROFILING == 1
                     , cl::vector<cl::Event*>* timings = nullptr
                     #endif // INOVESA_ENABLE_CLPROFILING
                     )
    {
        #if INOVESA_ENABLE_CLPROFILING == 1
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
    inline void
    enqueueWriteBuffer( const cl::Buffer& buffer
                      , cl_bool blocking
                      , cl::size_type src_offset
                      , cl::size_type size
                      , const void* ptr
                      , const cl::vector<cl::Event>* events = nullptr
                      , cl::Event* event = nullptr
                      #if INOVESA_ENABLE_CLPROFILING == 1
                      , cl::vector<cl::Event*>* timings = nullptr
                      #endif // INOVESA_ENABLE_CLPROFILING
                      )
    {
        #if INOVESA_ENABLE_CLPROFILING == 1
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

    inline void finish()
    {
        queue.finish();
    }

    inline void flush()
    {
        queue.flush();
    }

    #if INOVESA_ENABLE_CLPROFILING == 1
    void saveTimings(cl::vector<cl::Event *> *evts, std::string name);
    #endif // INOVESA_ENABLE_CLPROFILING

private:
    #if INOVESA_USE_CLFFT == 1
    clfftSetupData fft_setup;
    #endif // INOVESA_USE_CLFFT

    #if INOVESA_ENABLE_CLPROFILING == 1
    cl::vector<cl::Event*> timingsCopy;
    cl::vector<cl::Event*> timingsDFT;
    cl::vector<cl::Event*> timingsExecute;
    cl::vector<cl::Event*> timingsRead;
    cl::vector<cl::Event*> timingsWrite;
    #endif // INOVESA_ENABLE_CLPROFILING

    static const std::string custom_datatypes;

private: // helper functions
    static std::string datatype_aliases();

    static std::vector<cl_context_properties> properties(cl::Platform& platform,
                                                         bool glsharing);
};

#endif // INOVESA_USE_OPENCL
