// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

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

    ~Identity() noexcept override;

    /**
     * @brief apply copys data from in to out
     */
    void apply() override
    {
        #if INOVESA_USE_OPENCL == 1
        if (_oclh) {
            #if INOVESA_SYNC_CL == 1
            _in->syncCLMem(OCLH::clCopyDirection::cpu2dev);
            #endif // INOVESA_SYNC_CL
            _oclh->enqueueCopyBuffer( _in->data_buf, _out->data_buf
                                   , 0,0,sizeof(meshdata_t)*PhaseSpace::nxy
                                   #if INOVESA_ENABLE_CLPROFILING == 1
                                   , nullptr,nullptr
                                   , applySMEvents.get()
                                   # endif // INOVESA_ENABLE_CLPROFILING
                                   );
            _oclh->enqueueBarrier();
            #if INOVESA_SYNC_CL == 1
            _out->syncCLMem(OCLH::clCopyDirection::dev2cpu);
            #endif // INOVESA_SYNC_CL
        } else
        #endif // INOVESA_USE_OPENCL
        {
            meshdata_t* data_in = _in->getData();
            meshdata_t* data_out = _out->getData();
            std::copy_n(data_in,PhaseSpace::nb*PhaseSpace::nxy,data_out);
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
