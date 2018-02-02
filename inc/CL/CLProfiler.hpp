/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2017-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2017-2018: Karlsruhe Institute of Technology                 *
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

#ifndef CLPROFILER_HPP
#define CLPROFILER_HPP


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

#endif // CLPROFILER_HPP
