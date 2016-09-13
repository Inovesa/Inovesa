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

#ifndef IDENTITY_HPP
#define IDENTITY_HPP

#include "SM/SourceMap.hpp"

namespace vfps
{

class Identity : public SourceMap
{
public:
    Identity(PhaseSpace* in, PhaseSpace* out,
             const meshindex_t xsize, const meshindex_t ysize) :
        SourceMap(in, out, xsize, ysize, 0, 0) {}

    /**
     * @brief apply copys data from in to out
     */
    void apply()
    {
        #ifdef INOVESA_USE_CL
        if (OCLH::active) {
            #ifdef INOVESA_SYNC_CL
            _in->syncCLMem(clCopyDirection::cpu2dev);
            #endif // INOVESA_SYNC_CL
            OCLH::queue.enqueueCopyBuffer(_in->data_buf, _out->data_buf,
                                          0,0,sizeof(meshdata_t)*_size);
            #ifdef CL_VERSION_1_2
            OCLH::queue.enqueueBarrierWithWaitList();
            #else // CL_VERSION_1_2
            OCLH::queue.enqueueBarrier();
            #endif // CL_VERSION_1_2
            #ifdef INOVESA_SYNC_CL
            _out->syncCLMem(clCopyDirection::dev2cpu);
            #endif // INOVESA_SYNC_CL
        } else
        #endif // INOVESA_USE_CL
        {
            meshdata_t* data_in = _in->getData();
            meshdata_t* data_out = _out->getData();
            std::copy_n(data_in,_size,data_out);
        }
    }

    /**
     * @brief applyTo does nothing
     * @param x
     * @param y
     */
    PhaseSpace::position apply(const PhaseSpace::position pos) const
        { return pos; }
};

}

#endif // IDENTITY_HPP
