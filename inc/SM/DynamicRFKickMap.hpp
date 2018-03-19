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
    DynamicRFKickMap( std::shared_ptr<PhaseSpace> in
                    , std::shared_ptr<PhaseSpace> out
                    , const meshindex_t xsize
                    , const meshindex_t ysize
                    , const double dt
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
                    );

    ~DynamicRFKickMap() noexcept;

    /**
     * @brief apply updates RFKickMap before it is actually applied
     *
     * @todo Add OpenCL code path
     */
    void apply() override;

private:
    const meshaxis_t _addnoise;

    const meshaxis_t _mulnoise;

    const meshaxis_t _modampl;

    const meshaxis_t _modtimedelta;

    const uint32_t* _step;

    std::mt19937 _prng;
    std::normal_distribution<meshaxis_t> _dist;

    std::vector<std::array<meshaxis_t,2>> _modulation;

    /**
     * set up the KickMap with a (new) set of parameters
     */
    std::vector<std::array<meshaxis_t,2>> __calcModulation(uint32_t steps);

protected:
    /**
     * @brief update to current time step
     *
     * @todo For OpenCL, this should be done on device
     */
    void _update();
};

} // namespace vfps

#endif // DYNAMICRFKICKMAP_HPP
