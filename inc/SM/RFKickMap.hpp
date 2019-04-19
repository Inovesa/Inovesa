// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "SM/KickMap.hpp"

namespace vfps
{

class RFKickMap : public KickMap
{
public:
    RFKickMap(std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out
             , const meshaxis_t angle
             , const frequency_t f_RF
             , const InterpolationType it
             , const bool interpol_clamp
             , oclhptr_t oclh
             );


    RFKickMap(std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out
             , const timeaxis_t revolutionpart
             , const meshaxis_t V_RF
             , const frequency_t f_RF
             , const meshaxis_t V0
             , const InterpolationType it
             , const bool interpol_clamp
             , oclhptr_t oclh
             );

    ~RFKickMap() noexcept override;

protected:
    const bool _linear;

    const meshaxis_t _angle;

    const timeaxis_t _revolutionpart;

    const meshaxis_t _V_RF;

    const frequency_t _f_RF;

    const meshaxis_t _V0;

    const meshaxis_t _syncphase;

    const timeaxis_t _bl2phase;

protected:
    virtual void _calcKick( const meshaxis_t phase=0
                          , const meshaxis_t ampl=1);
};

} // namespace vfps

