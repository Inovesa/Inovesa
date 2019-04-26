// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include "SM/WakeKickMap.hpp"

vfps::WakeKickMap::WakeKickMap( std::shared_ptr<PhaseSpace> in
                              , std::shared_ptr<PhaseSpace> out
                              , const InterpolationType it
                              , const bool interpol_clamp
                              , oclhptr_t oclh
                              #if INOVESA_USE_OPENCL == 1 and INOVESA_USE_OPENGL == 1
                              , cl_GLuint glbuf
                              #endif // INOVESA_USE_OPENCL an INOVESA_USE_OPENGL
                              )  :
    KickMap( in,out,it,interpol_clamp,Axis::y,oclh)
    #if INOVESA_USE_OPENCL == 1 and INOVESA_USE_OPENGL == 1
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
