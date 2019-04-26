// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once
#if (INOVESA_USE_OPENCL == 1) && (INOVESA_ENABLE_CLPROFILING == 1)


#include <string>
#include <vector>


#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 110
#define CL_HPP_TARGET_OPENCL_VERSION 120
#include "local_cl.hpp"

namespace vfps
{

class CLTiming
{
public:
    CLTiming() = delete;

    CLTiming(const cl::Event& ev, std::string msg);

    const std::string msg;

    const cl_ulong submit;
    const cl_ulong queued;
    const cl_ulong start;
    const cl_ulong finish;
};

inline bool operator< (const CLTiming& lhs, const CLTiming& rhs)
    { return lhs.submit < rhs.submit; }

inline bool operator> (const CLTiming& lhs, const CLTiming& rhs)
    { return rhs < lhs; }

inline bool operator<=(const CLTiming& lhs, const CLTiming& rhs)
    { return !(lhs > rhs); }

inline bool operator>=(const CLTiming& lhs, const CLTiming& rhs)
    { return !(lhs < rhs); }

} // namespace vfps

#endif // INOVESA_USE_OPENCL, INOVESA_ENABLE_CLPROFILING
