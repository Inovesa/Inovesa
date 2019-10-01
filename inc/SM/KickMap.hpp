// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include "SM/SourceMap.hpp"

namespace vfps
{

/**
 * @brief The KickMap class allows to apply position dependent forces
 * and energy dependent displacements.
 *
 * For the KickMap the displacement is perpendicular to the gradient
 * on how large the displacement shall be.
 *
 */
class KickMap : public SourceMap
{
public:
    enum class Axis : bool {
        x=0, y=1
    };

public:
    KickMap( std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out
           , const meshindex_t xsize, const meshindex_t ysize
           , const InterpolationType it, const bool interpol_clamp
           , const Axis kd
           , oclhptr_t oclh
           );

    ~KickMap() noexcept;

public:
    const inline meshaxis_t* getForce() const
        { return _offset.data(); }

public:
    void apply() override;

    PhaseSpace::Position apply(PhaseSpace::Position pos) const override;

    #if INOVESA_USE_OPENCL == 1
    void syncCLMem(OCLH::clCopyDirection dir);
    #endif // INOVESA_USE_OPENCL

protected:
    /**
     * @brief _offset by one kick in units of mesh points
     */
    std::vector<meshaxis_t> _offset;

    #if INOVESA_USE_OPENCL == 1
    cl::Buffer _offset_clbuf;
    #endif // INOVESA_USE_OPENCL

    /**
     * @brief _kickdirection direction of the offset du to the kick
     */
    const Axis _kickdirection;

    /**
     * @brief _meshsize_kd size of the mesh in direction of the kick
     */
    #if INOVESA_USE_OPENCL == 1
    const cl_int _meshsize_kd;
    #else
    const meshindex_t _meshsize_kd;
    #endif // INOVESA_USE_OPENCL

    /**
     * @brief _meshsize_pd size of the mesh perpendicular to the kick
     */
    #if INOVESA_USE_OPENCL == 1
    const cl_int _meshsize_pd;
    #else
    const meshindex_t _meshsize_pd;
    #endif // INOVESA_USE_OPENCL

    /**
     * @brief updateSM
     *
     * does nothing when OpenCL is used
     */
    void updateSM();
};

}
