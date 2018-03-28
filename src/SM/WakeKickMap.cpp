/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
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

#include "SM/WakeKickMap.hpp"

vfps::WakeKickMap::WakeKickMap( std::shared_ptr<PhaseSpace> in
                              , std::shared_ptr<PhaseSpace> out
                              , const meshindex_t xsize
                              , const meshindex_t ysize
                              , const InterpolationType it
                              , const bool interpol_clamp
                              , oclhptr_t oclh
                              #if defined INOVESA_USE_OPENCL and defined INOVESA_USE_OPENGL
                              , cl_GLuint glbuf
                              #endif // INOVESA_USE_OPENCL an INOVESA_USE_OPENGL
                              )  :
    KickMap( in,out,xsize,ysize,it,interpol_clamp,Axis::y,oclh)
    #if defined INOVESA_USE_OPENCL and defined INOVESA_USE_OPENGL
    , _offset_glbuf(glbuf)
    #endif // INOVESA_USE_OPENCL and INOVESA_USE_OPENGL
{
}

vfps::WakeKickMap::~WakeKickMap() noexcept
#ifdef INOVESA_ENABLE_CLPROFILING
{
    saveTimings("WakeKickMap");
}
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING
