/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2015: Patrik Schönfeldt                                 *
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

#ifndef WAKEKICKMAP_HPP
#define WAKEKICKMAP_HPP

#include <array>
#include <fftw3.h>

#include "defines.hpp"
#include "IO/Display.hpp"
#include "HM/KickMap.hpp"
#include "ElectricField.hpp"
#include "Ruler.hpp"

using std::modf;

namespace vfps
{

/**
 * @brief The WakeKickMap class offers an option for one-dimensional kicks
 */
class WakeKickMap : public KickMap
{
public:
    /**
     * @brief WakeKickMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param wake wakefunction from -xsize to xsize-1, normalized in a way
     *             that plain multiplication with density gives an offset
     *
     * @todo interpolation when wake does not match PhaseSpace
     */
    WakeKickMap(PhaseSpace* in, PhaseSpace* out,
                const meshindex_t xsize, const meshindex_t ysize,
                const std::vector<std::pair<meshaxis_t,double>> wake,
                const InterpolationType it);

    /**
     * @brief WakeKickMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param csrimpedance
     * @param it
     */
    WakeKickMap(PhaseSpace* in, PhaseSpace* out,
                const meshindex_t xsize, const meshindex_t ysize,
                const vfps::ElectricField* csr,
                const InterpolationType it);

    ~WakeKickMap();

public:
    /**
     * @brief overloads KickMap::apply() to have a variable HeritageMap
     *
     * @todo currently uses phasespace in CPU Ram
     */
    void apply();

    inline const meshaxis_t* getWakeFunction() const
        { return _wakefunction; }

private:
    /**
     * @brief WakeKickMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param it
     */
    WakeKickMap(PhaseSpace* in, PhaseSpace* out,
                const meshindex_t xsize, const meshindex_t ysize,
                const InterpolationType it);

    /**
     * @brief _wakefunktion (normalized single particle) wake,
     *          sampled at 2*xsize positions [-xsize:+xsize]
     */
    meshaxis_t* const _wakefunction;

    const size_t _wakesize;
};

} // namespace VFPS

#endif // WAKEKICKMAP_HPP
