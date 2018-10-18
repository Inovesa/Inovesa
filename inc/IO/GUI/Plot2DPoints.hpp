/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2018: Patrik Schönfeldt                                      *
 * Copyright (c) 2018: Karlsruhe Institute of Technology                      *
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

#pragma once

#if INOVESA_USE_OPENGL == 1

#include "IO/GUI/Plot2DPrimitive.hpp"
#include "PS/PhaseSpace.hpp"

namespace vfps
{

class Plot2DPoints : public Plot2DPrimitive
{
public:
    Plot2DPoints() = delete;

    Plot2DPoints(std::array<float,3> rgb,
                 uint32_t nmeshcellsx,
                 uint32_t nmeshcellsy);

    virtual ~Plot2DPoints() noexcept = default;

    void update(const std::vector<PhaseSpace::Position>& points);

private:
    const uint32_t _nmeshcellsx;
    const uint32_t _nmeshcellsy;
};

} // namespace vfps

#endif // INOVESA_USE_OPENGL
