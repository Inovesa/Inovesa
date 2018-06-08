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

#include "SM/WakePotentialMap.hpp"

vfps::WakePotentialMap::WakePotentialMap( std::shared_ptr<PhaseSpace> in
                                        , std::shared_ptr<PhaseSpace> out
                                        , const vfps::meshindex_t xsize
                                        , const vfps::meshindex_t ysize
                                        , ElectricField* field
                                        , const InterpolationType it
                                        , bool interpol_clamp
                                        , oclhptr_t oclh
                                        )
  : WakeKickMap( in,out,xsize,ysize,it,interpol_clamp,oclh
               #if defined INOVESA_USE_OPENCL and defined INOVESA_USE_OPENGL
               , field->wakepotential_glbuf
               #endif // INOVESA_USE_OPENCL and INOVESA_USE_OPENGL
               )
  , _field(field)
{
}

vfps::WakePotentialMap::~WakePotentialMap() noexcept
#ifdef INOVESA_ENABLE_CLPROFILING
{
    saveTimings("WakePotentialMap");
}
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING

void vfps::WakePotentialMap::update()
{
    #ifdef INOVESA_USE_OPENCL
    if (_oclh) {
        _field->wakePotential();
        _oclh->enqueueCopyBuffer(_field->wakepotential_clbuf,_offset_clbuf,
                                0,0,sizeof(meshaxis_t)*_xsize);
        #ifdef INOVESA_SYNC_CL
        syncCLMem(OCLH::clCopyDirection::dev2cpu);
        #endif // INOVESA_SYNC_CL
    } else
    #endif // INOVESA_USE_OPENCL
    {
        std::copy_n(_field->wakePotential(),_xsize,_offset.data());
    }
    updateSM();
}
