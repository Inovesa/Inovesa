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

#ifndef WAKEFUNCTIONMAP_HPP
#define WAKEFUNCTIONMAP_HPP

#include "SM/WakeKickMap.hpp"

namespace vfps
{

class WakeFunctionMap : public WakeKickMap
{
public:
    /**
     * @brief WakeFunctionMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param fname
     */
    WakeFunctionMap(PhaseSpace* in, PhaseSpace* out,
                    const meshindex_t xsize, const meshindex_t ysize,
                    const std::string fname,
                    const double sigmaE, const double E0,
                    const double Ib, const double dt,
                    const InterpolationType it, const bool interpol_clamp);

    /**
     * @brief WakeFunctionMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param csrimpedance
     * @param it
     */
    WakeFunctionMap(PhaseSpace* in, PhaseSpace* out,
                    const meshindex_t xsize, const meshindex_t ysize,
                    const vfps::ElectricField* csr,
                    const InterpolationType it, const bool interpol_clamp);

    ~WakeFunctionMap();

public:
    inline const meshaxis_t* getWakeFunction() const
        { return _wakefunction; }

    /**
     * @brief update implements WakeKickMap
     *
     * @todo currently uses phasespace in CPU Ram
     */
    void update();

private:
    /**
     * @brief wakeFromFile
     * @param fname file name to read wake from
     * @return
     *
     * reads in a file and scales wake to internal units
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
     */
    WakeFunctionMap(PhaseSpace* in, PhaseSpace* out,
                    const meshindex_t xsize, const meshindex_t ysize,
                    const InterpolationType it, bool interpol_clamp);


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
