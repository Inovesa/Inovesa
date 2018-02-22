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

#ifndef PLOT2DLINE_HPP
#define PLOT2DLINE_HPP

#ifdef INOVESA_USE_OPENGL

#include <array>
#include <sstream>
#include <vector>

#include "IO/GUI/Plot2DPrimitive.hpp"

namespace vfps
{

class Plot2DLine : public Plot2DPrimitive
{
public:
    Plot2DLine() = delete;

    Plot2DLine(std::array<float,3> rgb);

    virtual ~Plot2DLine() noexcept = default;

    void update(const size_t npoints,
                const float* points,
                const bool vertical=false);


    void update(const size_t npoints,
                const double* points,
                const bool vertical=false);

private:
    float _max;
};

} // namespace vfps

#endif // INOVESA_USE_OPENGL

#endif // PLOT2DLINE_HPP
