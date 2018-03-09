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

#ifndef PLOT2DPRIMITIVE_HPP
#define PLOT2DPRIMITIVE_HPP

#ifdef INOVESA_USE_OPENGL


#include "IO/GUI/GUIElement.hpp"

#include <array>
#include <sstream>
#include <vector>

namespace vfps
{

class Plot2DPrimitive : public GUIElement
{
public:
    Plot2DPrimitive() = delete;

    Plot2DPrimitive(std::array<float,3> rgb,
                    size_t npoints,
                    GLenum primitivetype);

    virtual ~Plot2DPrimitive() noexcept;

    void draw() override;

    const cl::BufferGL& getCLBuffer() const
        { return _clBuffer; }

protected:
    std::vector<float> _points;

    size_t _npoints;

    GLuint vertexbuffer;
    GLuint position;

    cl::BufferGL _clBuffer;

private:
    const GLenum _primitivetype;

};

} // namespace vfps

#endif // INOVESA_USE_OPENGL

#endif // PLOT2DPRIMITIVE_HPP
