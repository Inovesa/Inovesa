// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#if INOVESA_USE_OPENGL == 1


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

protected:
    std::vector<float> _data;

    size_t _npoints;

    GLuint databuffer;
    GLuint dataID;

private:
    const GLenum _primitivetype;

};

} // namespace vfps

#endif // INOVESA_USE_OPENGL
