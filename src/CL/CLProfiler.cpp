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

#include "CL/CLProfiler.hpp"

vfps::CLTiming::CLTiming(const cl::Event& ev, std::string msg)
    : msg(msg)
    , submit(ev.getProfilingInfo<CL_PROFILING_COMMAND_SUBMIT>())
    , queued(ev.getProfilingInfo<CL_PROFILING_COMMAND_QUEUED>())
    , start(ev.getProfilingInfo<CL_PROFILING_COMMAND_START>())
    , finish(ev.getProfilingInfo<CL_PROFILING_COMMAND_END>())
{
}
