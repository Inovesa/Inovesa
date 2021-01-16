// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "SM/SourceMap.hpp"

/** \file
 *  \brief definitions of class KickMap
 */

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
    KickMap(std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out
           , const InterpolationType it, const bool interpol_clamp
           , const Axis kd
           , oclhptr_t oclh
           );

    ~KickMap() noexcept override;

public:
    /**
     * @brief getForce
     * @return pointer to beginning of _offset
     */
    const inline meshaxis_t* getForce() const
        { return _offset.data(); }

    /**
     * @brief swapOffset allows to replace _offset
     * @param offset
     */
    inline void swapOffset(std::vector<meshaxis_t>& offset)
        { _offset.swap(offset); updateSM(); }

public:
    void apply() override;

    void applyTo(PhaseSpace::Position& pos) const override;

    #if INOVESA_USE_OPENCL == 1
    void syncCLMem(OCLH::clCopyDirection dir);
    #endif // INOVESA_USE_OPENCL

protected:
    /**
     * @brief _offset by one kick in units of mesh points
     *
     * If kick is different for every bunch,
     * vector is used like a C-style ND array.
     *
     * @todo On the long run, a multi_array should be used.
     * If the kick is the same for every bunch,
     * the according dimentsion miht have just one entry.
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
     * @brief _lastbunch index of last bunch with individual KickMap
     */
    uint32_t _lastbunch;

    void updateSM();
};

}
