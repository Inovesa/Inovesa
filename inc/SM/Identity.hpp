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

#ifndef IDENTITY_HPP
#define IDENTITY_HPP

#include "SM/SourceMap.hpp"

namespace vfps
{

class Identity : public SourceMap
{
public:
    Identity( std::shared_ptr<PhaseSpace> in
            , std::shared_ptr<PhaseSpace> out
            , const meshindex_t xsize, const meshindex_t ysize
            , oclhptr_t oclh
            )
  : SourceMap( in, out, xsize, ysize, 0, 0, oclh)
{}

    ~Identity() noexcept;

    /**
     * @brief apply copys data from in to out
     */
    void apply() override
    {
        #ifdef INOVESA_USE_OPENCL
        if (_oclh) {
            #ifdef INOVESA_SYNC_CL
            _in->syncCLMem(OCLH::clCopyDirection::cpu2dev);
            #endif // INOVESA_SYNC_CL
            _oclh->enqueueCopyBuffer( _in->data_buf, _out->data_buf
                                   , 0,0,sizeof(meshdata_t)*_size
                                   #ifdef INOVESA_ENABLE_CLPROFILING
                                   , nullptr,nullptr
                                   , applySMEvents.get()
                                   # endif // INOVESA_ENABLE_CLPROFILING
                                   );
            _oclh->enqueueBarrier();
            #ifdef INOVESA_SYNC_CL
            _out->syncCLMem(OCLH::clCopyDirection::dev2cpu);
            #endif // INOVESA_SYNC_CL
        } else
        #endif // INOVESA_USE_OPENCL
        {
            meshdata_t* data_in = _in->getData();
            meshdata_t* data_out = _out->getData();
            std::copy_n(data_in,_size,data_out);
        }
    }

    /**
     * @brief applyTo does nothing
     */
    inline PhaseSpace::Position
    apply(const PhaseSpace::Position pos) const override
        { return pos; }
};

}

#endif // IDENTITY_HPP
