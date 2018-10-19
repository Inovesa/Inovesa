// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#if defined (INOVESA_USE_OPENCL) && defined (INOVESA_ENABLE_CLPROFILING)

#include "CL/CLProfiler.hpp"

vfps::CLTiming::CLTiming(const cl::Event& ev, std::string msg)
    : msg(msg)
    , submit(ev.getProfilingInfo<CL_PROFILING_COMMAND_SUBMIT>())
    , queued(ev.getProfilingInfo<CL_PROFILING_COMMAND_QUEUED>())
    , start(ev.getProfilingInfo<CL_PROFILING_COMMAND_START>())
    , finish(ev.getProfilingInfo<CL_PROFILING_COMMAND_END>())
{
}

#endif // INOVESA_USE_OPENCL, INOVESA_ENABLE_CLPROFILING
