// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
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
    WakeKickMap( std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out
               , const meshindex_t xsize, const meshindex_t ysize
               , const InterpolationType it, const bool interpol_clamp
               , oclhptr_t oclh
               #if defined INOVESA_USE_OPENCL and defined INOVESA_USE_OPENGL
               , vfps::clgluint glbuf
               #endif // INOVESA_USE_OPENCL and INOVESA_USE_OPENGL
               );

    ~WakeKickMap() noexcept;

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

