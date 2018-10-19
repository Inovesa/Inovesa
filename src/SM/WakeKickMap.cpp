// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

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
#if INOVESA_ENABLE_CLPROFILING == 1
{
    saveTimings("WakeKickMap");
}
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING
