// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#if INOVESA_USE_OPENGL == 1

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
    GUIElement() = delete;

    virtual ~GUIElement() noexcept;

    virtual void draw() =0;

    static uint_fast8_t glversion;

    inline bool getBufferShared() const
        { return buffer_shared; }

protected:
    GUIElement(bool buffer_shared);

    void compileShaders();

    std::string _fragmentshadercode;

    std::string _vertexshadercode;

    GLuint programID;

private:
    const bool buffer_shared;
};

} // namespace vfps

#endif // INOVESA_USE_OPENGL
