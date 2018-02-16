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

#include "SM/DriftMap.hpp"

#include <iostream>

vfps::DriftMap::DriftMap(std::shared_ptr<PhaseSpace> in,
                         std::shared_ptr<PhaseSpace> out,
                         const meshindex_t xsize,
                         const meshindex_t ysize,
                         const std::vector<meshaxis_t> slip,
                         const double E0,
                         const InterpolationType it,
                         const bool interpol_clamp)
    :
      KickMap(in,out,xsize,ysize,it,interpol_clamp,Axis::x)
{
    for(meshindex_t y=0; y<_ysize; y++) {
        _offset[y] = 0;
        for (size_t i=0; i<slip.size(); i++) {
            _offset[y] += slip[i]*_axis[1]->at(y)
                       *  std::pow(_axis[1]->at(y)*_axis[1]->scale()/E0,i);
        }
        _offset[y] /= _axis[1]->delta();
    }
    #ifdef INOVESA_USE_OPENCL
    if (OCLH::active) {
        syncCLMem(clCopyDirection::cpu2dev);
    }
    #endif // INOVESA_USE_OPENCL
    updateSM();
}

#ifdef INOVESA_ENABLE_CLPROFILING
vfps::DriftMap::~DriftMap() noexcept
{
    saveTimings("DrifMap");
}
#endif // INOVESA_ENABLE_CLPROFILING
