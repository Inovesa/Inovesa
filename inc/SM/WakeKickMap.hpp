// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include <array>
#include <fftw3.h>

#include "defines.hpp"
#include "IO/Display.hpp"
#include "SM/KickMap.hpp"
#include "PS/ElectricField.hpp"
#include "PS/Ruler.hpp"

namespace vfps
{

/**
 * @brief The WakeKickMap class offers an option for one-dimensional kicks
 */
class WakeKickMap : public KickMap
{
public:
    WakeKickMap(std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out
               , const InterpolationType it, const bool interpol_clamp
               , oclhptr_t oclh
               #if INOVESA_USE_OPENCL == 1 and INOVESA_USE_OPENGL == 1
               , vfps::clgluint glbuf
               #endif // INOVESA_USE_OPENCL and INOVESA_USE_OPENGL
               );

    ~WakeKickMap() noexcept override;

public:
    virtual void update()=0;

#if INOVESA_USE_OPENGL == 1
    vfps::clgluint getGLBuffer() const
        { return _offset_glbuf; }

protected:
    vfps::clgluint _offset_glbuf;
#endif // INOVESA_USE_OPENGL
};

} // namespace VFPS

