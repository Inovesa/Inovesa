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

#ifndef PLOT1DLINE_HPP
#define PLOT1DLINE_HPP

#ifdef INOVESA_USE_OPENGL

#include <array>
#include <sstream>
#include <vector>

#include "IO/GUI/GUIElement.hpp"

namespace vfps
{

class Plot1DLine : public GUIElement
{
public:
    enum class orientation : uint_fast16_t {
        horizontal,
        vertical
    };

public:
    Plot1DLine() = delete;

    Plot1DLine( std::array<float,3> rgb
              , size_t npoints
              , orientation orientation
              );

    virtual ~Plot1DLine() noexcept;

    void draw() override;

    void update(const float* points);


    void update(const double* points);

private:
    std::vector<float> _data;
    std::vector<float> _position;

    size_t _npoints;

    GLuint databuffer;
    GLuint dataID;

    GLuint posibuffer;
    GLuint posiID;

    void createPositionBuffer();

    float _max;

    orientation _orientation;
};

} // namespace vfps

#endif // INOVESA_USE_OPENGL

#endif // PLOT1DLINE_HPP
