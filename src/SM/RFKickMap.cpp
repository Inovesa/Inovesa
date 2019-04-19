// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include "SM/RFKickMap.hpp"

#include <boost/math/constants/constants.hpp>
using boost::math::constants::two_pi;

/**
 * @brief RFKickMap constructor for linear RF approximation
 */
vfps::RFKickMap::RFKickMap( std::shared_ptr<PhaseSpace> in
                          , std::shared_ptr<PhaseSpace> out
                          , const meshaxis_t angle
                          , const frequency_t f_RF
                          , const InterpolationType it
                          , const bool interpol_clamp
                          , oclhptr_t oclh
                          )
  : KickMap( in,out,it,interpol_clamp,Axis::y, oclh)
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

vfps::RFKickMap::RFKickMap( std::shared_ptr<PhaseSpace> in
                          , std::shared_ptr<PhaseSpace> out
                          , const timeaxis_t revolutionpart
                          , const meshaxis_t V_RF
                          , const frequency_t f_RF
                          , const meshaxis_t V0
                          , const InterpolationType it
                          , const bool interpol_clamp
                          , oclhptr_t oclh
                          )
  : KickMap( in,out,it,interpol_clamp,Axis::y, oclh)
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
