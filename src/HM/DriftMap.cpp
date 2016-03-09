/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Algorithms   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
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

#include "HM/DriftMap.hpp"

#include <iostream>

vfps::DriftMap::DriftMap(PhaseSpace *in, PhaseSpace *out,
                         const meshindex_t xsize,
                         const meshindex_t ysize,
                         const meshaxis_t angle,
                         const InterpolationType it)
    :
      KickMap(in,out,xsize,ysize,it,DirectionOfKick::x)
{
    for(meshindex_t y=0; y<_ysize; y++) {
        _offset[y] = std::tan(angle)*( static_cast<int>(y)
                                      -static_cast<int>(_ysize/2));
    }
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        syncCLMem(clCopyDirection::cpu2dev);
    }
    #endif // INOVESA_USE_CL
    updateHM();
}
