/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
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

vfps::WakePotentialMap::WakePotentialMap(vfps::PhaseSpace *in,
                                         vfps::PhaseSpace *out,
                                         const vfps::meshindex_t xsize,
                                         const vfps::meshindex_t ysize,
                                         ElectricField *field,
                                         const vfps::SourceMap::InterpolationType it,
                                         bool interpol_clamp):
    WakeKickMap(in,out,xsize,ysize,it,interpol_clamp),
    _field(field)
{
}

void vfps::WakePotentialMap::update()
{
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        _field->wakePotential();
        OCLH::queue.enqueueCopyBuffer(_field->_wakepotential_buf,_offset_buf,
                                      0,0,sizeof(meshaxis_t)*_xsize);
        #ifdef INOVESA_SYNC_CL
        syncCLMem(clCopyDirection::dev2cpu);
        #endif // INOVESA_SYNC_CL
    } else
    #endif // INOVESA_USE_CL
    {
        std::copy_n(_field->wakePotential(),_xsize,_offset.data());
    }
    updateHM();
}
