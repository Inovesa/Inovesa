// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Johannes Schestag
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include "SM/DynamicRFKickMap.hpp"

#include <cmath>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::two_pi;

vfps::DynamicRFKickMap::DynamicRFKickMap(std::shared_ptr<PhaseSpace> in
                                        , std::shared_ptr<PhaseSpace> out
                                        , const meshindex_t xsize
                                        , const meshindex_t ysize
                                        , const meshaxis_t angle
                                        , const double revolutionpart
                                        , const double f_RF
                                        , const meshaxis_t phasespread
                                        , const meshaxis_t amplspread
                                        , const meshaxis_t modampl
                                        , const double modtimeincrement
                                        , const uint32_t steps
                                        , const InterpolationType it
                                        , const bool interpol_clamp
                                        , oclhptr_t oclh
                                        )
  : RFKickMap( in,out,xsize,ysize,angle,f_RF,it,interpol_clamp, oclh)
  , _phasenoise(phasespread/std::sqrt(revolutionpart))
  , _amplnoise(amplspread/std::sqrt(revolutionpart))
  , _modampl(modampl)
  , _modtimedelta(two_pi<double>()*modtimeincrement)
  , _prng(std::mt19937(std::random_device{}()))
  , _dist(std::normal_distribution<meshaxis_t>(0, 1))
  , _next_modulation(__calcModulation(steps))
{
}

vfps::DynamicRFKickMap::DynamicRFKickMap( std::shared_ptr<PhaseSpace> in
                                        , std::shared_ptr<PhaseSpace> out
                                        , const meshindex_t xsize
                                        , const meshindex_t ysize
                                        , const double revolutionpart
                                        , const double V_RF
                                        , const double f_RF
                                        , const double V0
                                        , const meshaxis_t phasespread
                                        , const meshaxis_t amplspread
                                        , const meshaxis_t modampl
                                        , const double modtimeincrement
                                        , const uint32_t steps
                                        , const InterpolationType it
                                        , const bool interpol_clamp
                                        , oclhptr_t oclh
                                        )
  : RFKickMap( in,out,revolutionpart,V_RF,f_RF,V0,it,interpol_clamp,oclh )
  , _phasenoise(phasespread/std::sqrt(revolutionpart))
  , _amplnoise(amplspread/std::sqrt(revolutionpart))
  , _modampl(modampl)
  , _modtimedelta(two_pi<double>()*modtimeincrement)
  , _prng(std::mt19937(std::random_device{}()))
  , _dist(std::normal_distribution<meshaxis_t>(0, 1))
  , _next_modulation(__calcModulation(steps))
{
}

vfps::DynamicRFKickMap::~DynamicRFKickMap() noexcept
#if INOVESA_ENABLE_CLPROFILING == 1
{
    saveTimings("DynamicRFKickMap");
}
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING

std::queue<std::array<vfps::meshaxis_t,2>>
vfps::DynamicRFKickMap::__calcModulation(uint32_t steps)
{
    std::queue<std::array<meshaxis_t,2>> rv;

    for (uint32_t i=0; i<steps; i++) {
        meshaxis_t phasenoise = _dist(_prng)*_phasenoise;
        meshaxis_t amplnoise = _dist(_prng)*_amplnoise;

        meshaxis_t phasemod = _modampl*std::sin(_modtimedelta*(i));

        rv.emplace(std::array<meshaxis_t,2>{{ _syncphase+phasenoise+phasemod
                                            , 1+amplnoise}});
    }
    return rv;
}

void vfps::DynamicRFKickMap::_calcKick()
{
    RFKickMap::_calcKick( _next_modulation.front()[0]
                        , _next_modulation.front()[1]);
}

void vfps::DynamicRFKickMap::apply() {
    _calcKick();
    KickMap::apply();

    // move front entry from next to past
    _past_modulation.emplace_back(std::move(_next_modulation.front()));
    _next_modulation.pop();
}

std::vector<std::array<vfps::meshaxis_t,2>>
vfps::DynamicRFKickMap::getPastModulation()
{
    auto rv = std::move(_past_modulation);
    _past_modulation.clear();
    return rv;
}

