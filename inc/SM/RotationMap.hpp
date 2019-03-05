// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once


#include <array>

#include "defines.hpp"
#include "SM/SourceMap.hpp"
#include "IO/Display.hpp"

using std::modf;

namespace vfps
{

class RotationMap : public SourceMap
{
public:
    /**
     * @brief RotationMap
     * @param in data source
     * @param out data target
     * @param xsize horizontal size of thne mesh
     * @param ysize vertical size of thne mesh
     * @param angle to rotate the phase space
     * @param it number of points to use for (1D) interpolation
     * @param interpol_clamped turn clamping on and off
     * @param rt choose between real and Manhattan style
     * @param rotmapsize size of rotation map (can be 0, _size or _size/2)
     *
     * Saturation only makes sense with quadratic/cubic interpolation.
     * Enabling it with linear (or without)
     * interpolation currently crashes the program.
     */
    RotationMap(std::shared_ptr<PhaseSpace> in,
                std::shared_ptr<PhaseSpace> out,
                const meshindex_t xsize, const meshindex_t ysize,
                const meshaxis_t angle, const InterpolationType it,
                const bool interpol_clamped, const size_t rotmapsize,
                oclhptr_t oclh);

    #if INOVESA_ENABLE_CLPROFILING == 1
    ~RotationMap() noexcept;
    #else
    ~RotationMap() noexcept = default;
    #endif

    /**
     * @brief apply overrides HM::apply() by an optimized implementation
     */
    void apply() override;

    /**
     * @brief apply overrides HM::apply() by an optimized implementation
     *
     * @param pos position (in x/y coordinates)
     */
    void applyTo(PhaseSpace::Position& pos) const override;

private:
    const uint32_t _rotmapsize;

    const bool _clamp;

    void genHInfo(meshindex_t x0, meshindex_t y0, hi* myhinfo);

    const meshaxis_t _cos_dt;
    const meshaxis_t _sin_dt;

    #if INOVESA_USE_OPENCL == 1
    void genCode4SM4sat();

    void genCode4Rotation();

    cl_float2 rot;

    cl_int2 imgsize;

    cl_float2 zerobin;
    #endif // INOVESA_USE_OPENCL
};

}
