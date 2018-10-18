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
using boost::math::constants::two_pi;

vfps::RFKickMap::RFKickMap( std::shared_ptr<PhaseSpace> in
                          , std::shared_ptr<PhaseSpace> out
                          , const meshindex_t xsize
                          , const meshindex_t ysize
                          , const meshaxis_t angle
                          , const frequency_t f_RF
                          , const InterpolationType it
                          , const bool interpol_clamp
                          , oclhptr_t oclh
                          )
  : KickMap( in,out,xsize,ysize,it,interpol_clamp,Axis::y, oclh)
  , _linear(true)
  , _angle(angle)
  , _revolutionpart(0)
  , _V_RF(0)
  , _f_RF(f_RF)
  , _V0(0)
  , _syncphase(0)
  , _bl2phase(_axis[0]->scale("Meter")/physcons::c*_f_RF*two_pi<double>())
{
    _calcKick(_syncphase);
}

vfps::RFKickMap::RFKickMap(std::shared_ptr<PhaseSpace> in
                          , std::shared_ptr<PhaseSpace> out
                          , const meshindex_t xsize
                          , const meshindex_t ysize
                          , const timeaxis_t revolutionpart
                          , const meshaxis_t V_RF
                          , const frequency_t f_RF
                          , const meshaxis_t V0
                          , const InterpolationType it
                          , const bool interpol_clamp
                          , oclhptr_t oclh
                          )
  : KickMap( in,out,xsize,ysize,it,interpol_clamp,Axis::y, oclh)
  , _linear(false)
  , _angle(0)
  , _revolutionpart(revolutionpart)
  , _V_RF(V_RF)
  , _f_RF(f_RF)
  , _V0(V0)
  , _syncphase(std::asin(_V0/_V_RF))
  , _bl2phase(_axis[0]->scale("Meter")/physcons::c*_f_RF*two_pi<double>())
{
    _calcKick(_syncphase);
}

void vfps::RFKickMap::_calcKick(const meshaxis_t phase, const meshaxis_t ampl)
{
    if (_linear) {
        meshaxis_t phaseoffs = (_syncphase-phase);
        const meshaxis_t xcenter = _in->getAxis(0)->zerobin();
        for(meshindex_t x=0; x<_xsize; x++) {
            _offset[x] = std::tan(_angle)*(xcenter-x);
            _offset[x] += std::tan(_angle)*phaseoffs/_bl2phase/_axis[0]->delta();
            _offset[x] *= ampl;
        }
    } else {
        for(meshindex_t x=0; x<_xsize; x++) {
            _offset[x] = _revolutionpart*(-ampl*_V_RF
                       * std::sin(_axis[0]->at(x)*_bl2phase+phase)
                       + _V0)/ _axis[1]->delta()/_axis[1]->scale("ElectronVolt");
        }
    }

    #if INOVESA_USE_OPENCL == 1
    if (_oclh) {
        syncCLMem(OCLH::clCopyDirection::cpu2dev);
    }
    #endif // INOVESA_USE_OPENCL
    updateSM();
}

vfps::RFKickMap::~RFKickMap() noexcept
#if INOVESA_ENABLE_CLPROFILING == 1
{
    saveTimings("RFKickMap");
}
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING
