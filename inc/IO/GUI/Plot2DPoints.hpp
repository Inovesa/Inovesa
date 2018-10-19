// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

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
