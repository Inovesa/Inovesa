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

#ifndef WAKEFUNCTIONMAP_HPP
#define WAKEFUNCTIONMAP_HPP

#include "SM/WakeKickMap.hpp"

namespace vfps
{

/**
 * @brief The WakeFunctionMap class
 */
class WakeFunctionMap : public WakeKickMap
{
public:
    WakeFunctionMap() = delete;

    WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                   , std::shared_ptr<PhaseSpace> out
                   , const meshindex_t xsize, const meshindex_t ysize
                   , const std::string fname
                   , const double sigmaE, const double E0
                   , const double Ib, const double dt
                   , const InterpolationType it
                   , const bool interpol_clamp
                   , oclhptr_t oclh
                   );

    WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                   , std::shared_ptr<PhaseSpace> out
                   , const meshindex_t xsize, const meshindex_t ysize
                   , const ElectricField* csr
                   , const InterpolationType it
                   , const bool interpol_clamp
                   , oclhptr_t oclh
                   );

    ~WakeFunctionMap() noexcept;

public:
    inline const meshaxis_t* getWakeFunction() const
        { return _wakefunction; }

    /**
     * @brief update overrides WakeKickMap
     *
     * @todo currently uses phasespace in CPU Ram
     */
    void update() override;

private:
    /**
     * @brief _wakeFromFile reads in a file and scales wake to internal units
     * @param fname file name to read wake from
     */
    void _wakeFromFile(const std::string fname, const double scaling);


private:
    /**
     * @brief WakeFunctionMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param it
     * @param interpol_clamp
     * @param oclh
     *
     * @todo currently broken when used with CL/GL sharing
     */
    WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                   , std::shared_ptr<PhaseSpace> out
                   , const meshindex_t xsize, const meshindex_t ysize
                   , const InterpolationType it, bool interpol_clamp
                   , oclhptr_t oclh
                   );


    const Ruler<meshaxis_t> _xaxis;

    /**
     * @brief _wakefunktion (normalized single particle) wake,
     *          sampled at 2*xsize positions [-xsize:+xsize]
     */
    meshaxis_t* const _wakefunction;

    const size_t _wakesize;
};

} // namespace vfps

#endif // WAKEFUNCTIONMAP_HPP
