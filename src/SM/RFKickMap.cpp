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

#include "SM/RFKickMap.hpp"

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

vfps::RFKickMap::RFKickMap(std::shared_ptr<PhaseSpace> in
                          , std::shared_ptr<PhaseSpace> out
                          , const meshindex_t xsize
                          , const meshindex_t ysize
                          , const double dt
                          , const double V_RF
                          , const double f_RF
                          , const double V0
                          , const InterpolationType it
                          , const bool interpol_clamp
                          , oclhptr_t oclh
                          )
  : KickMap( in,out,xsize,ysize,it,interpol_clamp,Axis::y, oclh)
{
    auto syncphase = std::asin(V0/V_RF);
    auto bl2phase = _axis[0]->scale()/physcons::c*f_RF/pi<double>();

    for(meshindex_t x=0; x<_xsize; x++) {
        _offset[x] = dt*(-V_RF*std::sin(_axis[0]->at(x)*bl2phase+syncphase)+V0);
        _offset[x] /= _axis[0]->delta()*_axis[0]->scale();
    }
    #ifdef INOVESA_USE_OPENCL
    if (_oclh) {
        syncCLMem(clCopyDirection::cpu2dev);
    }
    #endif // INOVESA_USE_OPENCL
    updateSM();
}

vfps::RFKickMap::~RFKickMap() noexcept
#ifdef INOVESA_ENABLE_CLPROFILING
{
    saveTimings("RFKickMap");
}
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING
