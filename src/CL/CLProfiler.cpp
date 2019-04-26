// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#if (INOVESA_USE_OPENCL == 1) && (INOVESA_ENABLE_CLPROFILING == 1)

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
