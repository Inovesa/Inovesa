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

#ifndef DYNAMICRFKICKMAP_HPP
#define DYNAMICRFKICKMAP_HPP

#include "SM/RFKickMap.hpp"

#include <random>

namespace vfps
{

/*!
 * @brief Provides the RFKickMap with a random zeroth and first order contribution
 */

class DynamicRFKickMap : public RFKickMap
{
public:
    /*!
     * @param angle     the RFKickMap rotation angle
     * @param s_phase   width of the additive instability in y units
     * @param s_peak    width of the relative multiplicative instability
     */
    DynamicRFKickMap(std::shared_ptr<PhaseSpace> in,
                     std::shared_ptr<PhaseSpace> out,
                     const meshindex_t xsize, const meshindex_t ysize,
                     const meshaxis_t angle,
                     const meshaxis_t add,
                     const meshaxis_t multiply,
                     const InterpolationType it,
                     const bool interpol_clamp);

    ~DynamicRFKickMap() noexcept;

    /**
     * @brief apply updates RFKickMap before it is actually applied
     *
     * @todo Add OpenCL code path
     */
    void apply() override;

private:
    const meshaxis_t _add, _multiply;
    std::mt19937 _prng;
    std::normal_distribution<meshaxis_t> _dist;

    std::vector<meshaxis_t> _mean;

    /// set up the KickMap with a new set of random parameters
    void reset();
};

} // namespace vfps

#endif // DYNAMICRFKICKMAP_HPP
