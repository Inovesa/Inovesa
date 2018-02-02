#include "CL/CLProfiler.hpp"

vfps::CLTiming::CLTiming(const cl::Event& ev, std::string msg)
    : msg(msg)
    , submit(ev.getProfilingInfo<CL_PROFILING_COMMAND_SUBMIT>())
    , queued(ev.getProfilingInfo<CL_PROFILING_COMMAND_QUEUED>())
    , start(ev.getProfilingInfo<CL_PROFILING_COMMAND_START>())
    , finish(ev.getProfilingInfo<CL_PROFILING_COMMAND_END>())
{
}
