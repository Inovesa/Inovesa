/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2017: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2017: Karlsruhe Institute of Technology                 *
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

#ifndef GUIELEMENT_HPP
#define GUIELEMENT_HPP

#ifdef INOVESA_USE_GUI

//forward declaration
namespace vfps {
class GUIElement;
}

#include <exception>
#include <fstream>
#include <string>
#include <vector>

#include <GL/glew.h>
#if GLFW_VERSION_MAJOR == 3
#include <GLFW/glfw3.h>
#else
#include <GL/glfw.h>
#endif

#include "PS/PhaseSpace.hpp"

namespace vfps
{

class GUIElement
{
public:
    GUIElement();

    virtual ~GUIElement() noexcept;

    virtual void draw() =0;

    static uint_fast8_t glversion;

protected:
    void compileShaders();

    std::string _fragmentshadercode;

    std::string _vertexshadercode;

    GLuint programID;
};

} // namespace vfps

#endif // INOVESA_USE_GUI

#endif // GUIELEMENT_HPP
