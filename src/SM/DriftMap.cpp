// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#include "SM/DriftMap.hpp"

#include <iostream>

vfps::DriftMap::DriftMap( std::shared_ptr<PhaseSpace> in
                        , std::shared_ptr<PhaseSpace> out
                        , const meshindex_t xsize
                        , const meshindex_t ysize
                        , const std::vector<meshaxis_t>& slip
                        , const double E0
                        , const InterpolationType it
                        , const bool interpol_clamp
                        , oclhptr_t oclh
                        )
  : KickMap( in,out,xsize,ysize,it,interpol_clamp,Axis::x,oclh)
{
    for(meshindex_t y=0; y<_ysize; y++) {
        _offset[y] = 0;
        for (size_t i=0; i<slip.size(); i++) {
            _offset[y] += slip[i]*_axis[1]->at(y)
                       *  std::pow(_axis[1]->at(y)*_axis[1]->scale("ElectronVolt")/E0,i);
        }
        _offset[y] /= _axis[0]->delta();
    }

    #if INOVESA_USE_OPENCL == 1
    if (_oclh) {
        syncCLMem(OCLH::clCopyDirection::cpu2dev);
    }
    #endif // INOVESA_USE_OPENCL
    updateSM();
}

#if INOVESA_ENABLE_CLPROFILING == 1
vfps::DriftMap::~DriftMap() noexcept
{
    saveTimings("DrifMap");
}
#endif // INOVESA_ENABLE_CLPROFILING
