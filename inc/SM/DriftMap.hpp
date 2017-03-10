/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2017: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2017: Karlsruhe Institute of Technology                 *
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

#ifndef DRIFTMAP_HPP
#define DRIFTMAP_HPP

#include "SM/KickMap.hpp"

namespace vfps
{

class DriftMap : public KickMap
{
public:
    DriftMap(std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out,
             const meshindex_t xsize, const meshindex_t ysize,
             const std::vector<meshaxis_t> slip, const double E0,
             const InterpolationType it, const bool interpol_clamp);

    ~DriftMap()
    #ifdef INOVESA_ENABLE_CLPROFILING
        { std::cout << "~DriftMap() -> "; }
    #else
    = defaut;
    #endif // INOVESA_ENABLE_CLPROFILING
};

} // namespace fvps

#endif // DRIFTMAP_HPP
