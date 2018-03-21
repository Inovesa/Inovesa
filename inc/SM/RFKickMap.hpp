/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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

#ifndef RFKICKMAP_HPP
#define RFKICKMAP_HPP

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

#endif // RFKICKMAP_HPP
