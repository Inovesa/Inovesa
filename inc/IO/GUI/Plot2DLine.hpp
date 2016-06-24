/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
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

#ifdef INOVESA_USE_GUI

#include <array>
#include <sstream>
#include <vector>

#include "IO/GUI/GUIElement.hpp"

namespace vfps
{

class Plot2DLine : public GUIElement
{
public:
    Plot2DLine(std::array<float,3> rgb);

    ~Plot2DLine();

    void updateLine(const size_t npoints, const float* points);

    void draw();

private:
    size_t _npoints;

    std::vector<float> _line;

    GLuint vertexbuffer;
    GLuint position;
};

} // namespace vfps

#endif // INOVESA_USE_GUI

#endif // PLOT2DLINE_HPP
