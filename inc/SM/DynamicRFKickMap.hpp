// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Johannes Schestag
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "SM/RFKickMap.hpp"

#include <queue>
#include <random>

namespace vfps
{

/*!
 * @brief Provides the RFKickMap with a random zeroth and first order contribution
 */

class DynamicRFKickMap : public RFKickMap
{
public:
    /**
     * @brief constructor for linear RF
     *
     * @param angle the RFKickMap rotation angle
     * @param phasespread width of the additive instability in rad
     */
    DynamicRFKickMap( std::shared_ptr<PhaseSpace> in
                    , std::shared_ptr<PhaseSpace> out
                    , const meshindex_t xsize
                    , const meshindex_t ysize
                    , const meshaxis_t angle
                    , const double revolutionpart
                    , const double f_RF
                    , const meshaxis_t phasespread
                    , const meshaxis_t mulnoise
                    , const meshaxis_t modampl
                    , const double modtimeincrement
                    , const uint32_t steps
                    , const InterpolationType it
                    , const bool interpol_clamp
                    , oclhptr_t oclh
                    );
    /**
     * @brief constructor for sinusoidal RF
     *
     * @param phasespread width of the additive instability in rad
     */
    DynamicRFKickMap(std::shared_ptr<PhaseSpace> in
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
                    );

    ~DynamicRFKickMap() noexcept override;

    /**
     * @brief apply updates RFKickMap before it is actually applied
     *
     * @todo Add OpenCL code path
     */
    void apply() override;

    /**
     * @brief getPastModulation
     * @return modulations (phase,amplitude) since last call
     */
    std::vector<std::array<meshaxis_t,2>> getPastModulation();

private:
    const meshaxis_t _phasenoise;

    const meshaxis_t _amplnoise;

    const meshaxis_t _modampl;

    const meshaxis_t _modtimedelta;

    std::mt19937 _prng;
    std::normal_distribution<meshaxis_t> _dist;

    std::queue<std::array<meshaxis_t,2>> _next_modulation;

    std::vector<std::array<meshaxis_t,2>> _past_modulation;

    /**
     * set up the KickMap with a (new) set of parameters
     */
    std::queue<std::array<meshaxis_t,2>> __calcModulation(uint32_t steps);

protected:
    /**
     * @brief update to current time step
     *
     * @todo For OpenCL, this should be done on device
     */
    void _calcKick();
};

} // namespace vfps

