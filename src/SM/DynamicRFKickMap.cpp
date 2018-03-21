/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
 * Copyright (c) 2017-2018: Johannes Schestag                                 *
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

#include "SM/DynamicRFKickMap.hpp"

#include <cmath>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::two_pi;

vfps::DynamicRFKickMap::DynamicRFKickMap( std::shared_ptr<PhaseSpace> in
                                        , std::shared_ptr<PhaseSpace> out
                                        , const meshindex_t xsize
                                        , const meshindex_t ysize
                                        , const meshaxis_t angle
                                        , const double f_RF
                                        , const meshaxis_t addnoise
                                        , const meshaxis_t mulnoise
                                        , const meshaxis_t modampl
                                        , const double modtimeincrement
                                        , const uint32_t* step
                                        , const uint32_t steps
                                        , const InterpolationType it
                                        , const bool interpol_clamp
                                        , oclhptr_t oclh
                                        )
    : RFKickMap( in,out,xsize,ysize,angle,f_RF,it,interpol_clamp, oclh)
    , _addnoise(addnoise)
    , _mulnoise(mulnoise)
    , _modampl(modampl)
    , _modtimedelta(two_pi<double>()*modtimeincrement)
    , _step(step)
    , _prng(std::mt19937(std::random_device{}()))
    , _dist(std::normal_distribution<meshaxis_t>(0, 1))
    , _modulation(__calcModulation(steps))
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
                                        , const meshaxis_t addnoise
                                        , const meshaxis_t mulnoise
                                        , const meshaxis_t modampl
                                        , const double modtimeincrement
                                        , const uint32_t* step
                                        , const uint32_t steps
                                        , const InterpolationType it
                                        , const bool interpol_clamp
                                        , oclhptr_t oclh
                                        )
  : RFKickMap( in,out,xsize,ysize
             , revolutionpart,V_RF,f_RF,V0,it,interpol_clamp,oclh )
  , _addnoise(addnoise)
  , _mulnoise(mulnoise)
  , _modampl(modampl)
  , _modtimedelta(two_pi<double>()*modtimeincrement)
  , _step(step)
  , _prng(std::mt19937(std::random_device{}()))
  , _dist(std::normal_distribution<meshaxis_t>(0, 1))
  , _modulation(__calcModulation(steps))
{
}

vfps::DynamicRFKickMap::~DynamicRFKickMap() noexcept
#ifdef INOVESA_ENABLE_CLPROFILING
{
    saveTimings("DynamicRFKickMap");
}
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING

std::vector<std::array<vfps::meshaxis_t,2>>
vfps::DynamicRFKickMap::__calcModulation(uint32_t steps)
{
    std::vector<std::array<meshaxis_t,2>> rv;
    rv.reserve(steps);

    for (uint32_t i=0; i<steps; i++) {
        meshaxis_t phasenoise = _dist(_prng)*_addnoise;
        meshaxis_t amplnoise = _dist(_prng)*_mulnoise;

        meshaxis_t phasemod = _modampl*std::sin(_modtimedelta*(i));

        rv.push_back({{ _syncphase+phasenoise+phasemod,1+amplnoise}});
    }
    return rv;
}

void vfps::DynamicRFKickMap::_calcKick()
{
    RFKickMap::_calcKick(_modulation[*_step][0],_modulation[*_step][1]);
}

void vfps::DynamicRFKickMap::apply() {
    _calcKick();
    KickMap::apply();
}


