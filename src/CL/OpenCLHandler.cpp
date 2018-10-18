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

#if INOVESA_USE_OPENCL == 1

#include "CL/OpenCLHandler.hpp"
#include "IO/Display.hpp"
#ifdef __linux__
#include <GL/glx.h>
#endif
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl_gl.h>
#include <OpenCL/cl_gl_ext.h>
#include <OpenGL/CGLDevice.h>
#include <OpenGL/CGLCurrent.h>
#else
#include <CL/cl_gl.h>
#endif

OCLH::OCLH( uint32_t device, bool glsharing)
  : ogl_sharing(glsharing)
{
    cl::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    uint32_t devicescount = 0;

    for (auto p : platforms) {
        try {
            context = cl::Context( CL_DEVICE_TYPE_ALL
                                 , properties(p,ogl_sharing).data());
        } catch (cl::Error& e) {
            ogl_sharing = false;
            if (e.err() == CL_INVALID_PROPERTY) {
                context = cl::Context( CL_DEVICE_TYPE_ALL
                                     , properties(p,ogl_sharing).data());
            } else {
                throw e;
            }
        }
        _devices = context.getInfo<CL_CONTEXT_DEVICES>();
        if (devicescount + _devices.size() <= device) {
            // device is on later platform
            devicescount += _devices.size();
        } else {
            _platform = p;
            // subtract devices on previous platforms
            _device = _devices[device-devicescount];
            break;
        }
    }

    #if INOVESA_ENABLE_CLPROFILING == 1
    queue = cl::CommandQueue(context,_device,CL_QUEUE_PROFILING_ENABLE);
    #else
    queue = cl::CommandQueue(context,_device);
    #endif // INOVESA_ENABLE_CLPROFILING

    devicetype = _device.getInfo<CL_DEVICE_TYPE>();

    #if INOVESA_USE_OPENGL == 1
    // cl_VENDOR_gl_sharing is present, when string contains the substring
    if (_device.getInfo<CL_DEVICE_EXTENSIONS>().find("_gl_sharing")
            == std::string::npos) {
        ogl_sharing = false;
    }
    #endif // INOVESA_USE_OPENGL

    #if INOVESA_USE_CLFFT == 1
    clfftInitSetupData(&fft_setup);
    clfftSetup(&fft_setup);
    #endif // INOVESA_USE_CLFFT

    #if INOVESA_ENABLE_CLPROFILING == 1
    // place initial marker
    queue.enqueueMarkerWithWaitList(nullptr,&init);
    #endif // INOVESA_ENABLE_CLPROFILING


    vfps::Display::printText("Initialized \""
                             + _device.getInfo<CL_DEVICE_NAME>()
                             + "\" (on platform \""
                             + _platform.getInfo<CL_PLATFORM_NAME>()
                             + "\") for use with OpenCL.");
    if (ogl_sharing) {
        vfps::Display::printText("Sharing between OpenCL "
                                 "and OpenGL is active.");
    }
}

cl::Program OCLH::prepareCLProg(std::string code)
{
    code = datatype_aliases()+custom_datatypes+code;
    cl::vector<std::string> codevec;
    codevec.push_back(code);
    cl::Program::Sources source(codevec);
    cl::Program p(context, source);
    try {
        // empty for compatibility reasons.
        std::string OCLBuildOpts("");

        p.build(_devices,OCLBuildOpts.c_str());
    } catch (cl::Error &e) {
        std::cerr << e.what() << std::endl;
        std::cout << "===== OpenCL Code =====\n"
                                << code << std::endl;
    #if DEBUG == 1
        throw e;
    }
    #endif
        std::cout << "===== OpenCL Build Log =====\n"
                  << p.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_device)
                  << std::endl;
    #ifndef DEBUG
        throw e;
    }
    #endif

return p;
}

#if INOVESA_ENABLE_CLPROFILING == 1
void OCLH::saveProfilingInfo(std::string fname)
{
    saveTimings(&timingsCopy,"MiscCopy");
    saveTimings(&timingsDFT,"DFTrafo");
    saveTimings(&timingsExecute,"MiscExec");
    saveTimings(&timingsRead,"MiscRead");
    saveTimings(&timingsWrite,"MiscWrite");

    std::ofstream timefile(fname);
    cl_ulong starttime(init.getProfilingInfo<CL_PROFILING_COMMAND_SUBMIT>());
    timingInfo.sort();
    timefile << "submit"
     << '\t' << "queued"
     << '\t' << "start"
     << '\t' << "finish"
     << '\t' << "type"
     << std::endl;
    for (auto ev : timingInfo) {
        timefile << ev.submit-starttime
         << '\t' << ev.queued-starttime
         << '\t' << ev.start-starttime
         << '\t' << ev.finish-starttime
         << '\t' << ev.msg
         << std::endl;
    }
}
#endif // INOVESA_ENABLE_CLPROFILING

#if INOVESA_ENABLE_CLPROFILING == 1
void OCLH::saveTimings(cl::vector<cl::Event*>* evts, std::string name)
{
    queue.finish();
    for (auto ev : *evts) {
        try {
            OCLH::timingInfo.push_back(vfps::CLTiming(*ev,name));
        } catch (...) {
            std::cerr << "Error in " << name << std::endl;
        }
    }
}
#endif // INOVESA_ENABLE_CLPROFILING

OCLH::~OCLH()
{
    #if INOVESA_ENABLE_CLPROFILING == 1
    saveProfilingInfo("inovesa-timings.txt");
    #else
    OCLH::queue.finish();
    #endif

    #if INOVESA_USE_CLFFT == 1
    clfftTeardown();
    #endif // INOVESA_USE_CLFFT

    #if INOVESA_ENABLE_CLPROFILING == 1
    for (auto ev : timingsCopy) {
        delete ev;
    }
    for (auto ev : timingsDFT) {
        delete ev;
    }
    for (auto ev : timingsExecute) {
        delete ev;
    }
    for (auto ev : timingsRead) {
        delete ev;
    }
    for (auto ev : timingsWrite) {
        delete ev;
    }
    #endif // INOVESA_ENABLE_CLPROFILING
}

void OCLH::listCLDevices()
{
    cl::vector<cl::Platform> platforms;
    std::cout << "OpenCL device options available on this computer:" << std::endl
              << " 0: (Do not use OpenCL.)" << std::endl;
    try {
        cl::Platform::get(&platforms);
    } catch (cl::Error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return;
    }

    uint32_t devicescount = 0;

    for (auto p : platforms) {
        std::string available_clversion;
        p.getInfo(CL_PLATFORM_VERSION,&available_clversion);
        std::string platformname;
        p.getInfo(CL_PLATFORM_NAME,&platformname);
        std::cout << "On platform " << platformname
                  << " (" << available_clversion << ")" << ":" << std::endl;
        cl::Context tmp_context = cl::Context(CL_DEVICE_TYPE_ALL, properties(p,false).data());

        cl::vector<cl::Device> tmp_devices = tmp_context.getInfo<CL_CONTEXT_DEVICES>();
        for (unsigned int d=0; d<tmp_devices.size(); d++) {
            devicescount++;
            std::string tmpdevicetype;
            switch(tmp_devices[d].getInfo<CL_DEVICE_TYPE>()) {
            case CL_DEVICE_TYPE_CPU:
                tmpdevicetype = "CPU";
                break;
            case CL_DEVICE_TYPE_GPU:
                tmpdevicetype = "GPU";
                break;
            default:
                tmpdevicetype = "unknown";
                break;
            }
            std::cout << " " << devicescount << ": "
                      << tmpdevicetype << ", "
                      << tmp_devices[d].getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>()/0x100000 << " MiB, "
                      << tmp_devices[d].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << "CU "
                      << "(" << tmp_devices[d].getInfo<CL_DEVICE_NAME>() << ")"
                      #if DEBUG == 1
                      << std::endl
                      << "\ton \"" <<  p.getInfo<CL_PLATFORM_NAME>() << "\""
                      << " by \"" <<  p.getInfo<CL_PLATFORM_VENDOR>() << '\"'
                      << std::endl
                      << "\toffering \"" << tmp_devices[d].getInfo<CL_DEVICE_VERSION>() << '\"'
                      << " with " << tmp_devices[d].getInfo<CL_DEVICE_EXTENSIONS>()
                      #endif // DEBUG
                      << std::endl;
        }
    }
}

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
    } else if (std::is_same<vfps::meshdata_t,double>::value) {
    code +=
        "typedef double data_t;\n"
        "typedef double2 data2_t;\n"
        "typedef double3 data3_t;\n"
        "typedef double4 data4_t;\n"
        "double mult(double x, double y);"
        "double mult(double x, double y) { return x*y; }\n";
    }
    return code;
}

std::vector<cl_context_properties> OCLH::properties( cl::Platform& platform
                                                   , bool glsharing)
{
    std::vector<cl_context_properties> rv;

    #if INOVESA_USE_OPENGL == 1
    if (glsharing) {
        #ifdef __linux__
        rv = {
            CL_GL_CONTEXT_KHR,
            reinterpret_cast<cl_context_properties>(glXGetCurrentContext()),
            CL_GLX_DISPLAY_KHR,
            reinterpret_cast<cl_context_properties>(glXGetCurrentDisplay()),
            CL_CONTEXT_PLATFORM,
            reinterpret_cast<cl_context_properties>(platform()),
            0
        };
        #elif defined _WIN32
        rv = {
            CL_GL_CONTEXT_KHR,
            reinterpret_cast<cl_context_properties>(wglGetCurrentContext()),
            CL_WGL_HDC_KHR,
            reinterpret_cast<cl_context_properties>(wglGetCurrentDC()),
            CL_CONTEXT_PLATFORM,
            reinterpret_cast<cl_context_properties>(platform()),
            0
        };
        #elif defined TARGET_OS_MAC
        CGLContextObj glContext = CGLGetCurrentContext();
        CGLShareGroupObj shareGroup = CGLGetShareGroup(glContext);
        rv = {
            CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE,
            reinterpret_cast<cl_context_properties>(shareGroup),
            0
        };
        #endif
    } else
    #endif // INOVESA_USE_OPENGL
    {
        rv = {
            CL_CONTEXT_PLATFORM,
            reinterpret_cast<cl_context_properties>(platform()),
            0
        };
    }

    return rv;
}

#endif // INOVESA_USE_OPENCL
