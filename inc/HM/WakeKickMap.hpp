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

namespace vfps
{

/**
 * @brief The WakeKickMap class offers an option for one-dimensional kicks
 */
class WakeKickMap : public KickMap
{
public:
    WakeKickMap(vfps::PhaseSpace* in, vfps::PhaseSpace* out,
                const meshindex_t xsize, const meshindex_t ysize,
                const InterpolationType it);

    ~WakeKickMap();

public:
    /**
     * @brief overloads KickMap::apply() to have a variable HeritageMap
     *
     * @todo currently uses phasespace in CPU Ram
     */
    void apply();

    virtual void update()=0;
};

} // namespace VFPS

#endif // WAKEKICKMAP_HPP
