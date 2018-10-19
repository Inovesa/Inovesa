// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include "SM/KickMap.hpp"

namespace vfps
{

class RFKickMap : public KickMap
{
public:
    /**
     * @brief RFKickMap constructor for linear RF approximation
     */
    RFKickMap( std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out
             , const meshindex_t xsize, const meshindex_t ysize
             , const meshaxis_t angle
             , const frequency_t f_RF
             , const InterpolationType it
             , const bool interpol_clamp
             , oclhptr_t oclh
             );


    RFKickMap( std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out
             , const meshindex_t xsize, const meshindex_t ysize
             , const timeaxis_t revolutionpart
             , const meshaxis_t V_RF
             , const frequency_t f_RF
             , const meshaxis_t V0
             , const InterpolationType it
             , const bool interpol_clamp
             , oclhptr_t oclh
             );

    ~RFKickMap() noexcept;

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

