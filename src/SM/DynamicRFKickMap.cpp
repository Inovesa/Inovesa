/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2016: Karlsruhe Institute of Technology                 *
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


vfps::DynamicRFKickMap::DynamicRFKickMap(std::shared_ptr<PhaseSpace> in,
                           std::shared_ptr<PhaseSpace> out,
                           const meshindex_t xsize,
                           const meshindex_t ysize,
                           const meshaxis_t angle,
                           const meshaxis_t add,
                           const meshaxis_t multiply,
                           const InterpolationType it,
                           const bool interpol_clamp)
    :
      RFKickMap(in,out,xsize,ysize,angle,it,interpol_clamp),
      _add(add),
      _multiply(multiply),
      _prng(std::mt19937(std::random_device{}())),
      _dist(std::normal_distribution<meshaxis_t>(0, 1)),
      _mean(_offset)
{
}

vfps::DynamicRFKickMap::~DynamicRFKickMap()
#ifdef INOVESA_ENABLE_CLPROFILING
    { std::cout << "~DynamicRFKickMap() -> "; }
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING

void vfps::DynamicRFKickMap::reset() {
    meshaxis_t mul = _dist(_prng)*_multiply;
    meshaxis_t add = _dist(_prng)*_add;

    _offset = _mean;

    std::for_each(_offset.begin(), _offset.end(),
            [&mul, &add](auto &it) { it = fma(it, 1 + mul, add); }
            );

    #ifdef INOVESA_USE_OPENCL
    if (OCLH::active) {
        syncCLMem(clCopyDirection::cpu2dev);
    } else
    #endif // INOVESA_USE_OPENCL
    {
        updateSM();
    }
}

void vfps::DynamicRFKickMap::apply() {
    reset();
    KickMap::apply();
}


