// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#if INOVESA_USE_OPENGL == 1

#include <array>
#include <sstream>
#include <vector>

#include "IO/GUI/GUIElement.hpp"

namespace vfps
{

class Plot1DLine : public GUIElement
{
public:
    enum class Orientation : uint_fast16_t {
        horizontal,
        vertical
    };

public:
    Plot1DLine() = delete;

    /**
     * @brief Plot1DLine
     * @param rgb
     * @param npoints
     * @param orientation
     * @param databuffer takes ownership or creates (if 0)
     */
    Plot1DLine( std::array<float,3> rgb
              , size_t npoints
              , Orientation orientation
              , GLuint databuffer
              );

    virtual ~Plot1DLine() noexcept;

    void draw() override;

    /**
     * @brief update data for line to be renedered (meant for non-OpenCL mode)
     * @param points array of _npoints to be rendered as line
     *
     * not needed when shared buffer is updated from another location
     */
    void update(const float* points);

private:
    std::vector<float> _data;
    std::vector<float> _position;

    size_t _npoints;

    GLuint _databuffer;
    GLuint dataID;

    GLuint _posibuffer;
    GLuint posiID;

    void createPositionBuffer();

    Orientation _orientation;
};

} // namespace vfps

#endif // INOVESA_USE_OPENGL
