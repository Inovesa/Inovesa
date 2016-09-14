/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
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

#ifndef KICKMAP_HPP
#define KICKMAP_HPP

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
    KickMap(PhaseSpace* in, PhaseSpace* out,
            const meshindex_t xsize, const meshindex_t ysize,
            const InterpolationType it, const bool interpol_clamp,
            const Axis kd);

    ~KickMap();

public:
    const inline meshaxis_t* getForce() const
        { return _offset.data(); }

public:
    void apply();

    PhaseSpace::Position apply(PhaseSpace::Position pos) const;

    #ifdef INOVESA_USE_CL
    void syncCLMem(clCopyDirection dir);
    #endif // INOVESA_USE_CL

protected:
    /**
     * @brief _offset by one kick in units of mesh points
     */
    std::vector<meshaxis_t> _offset;

    #ifdef INOVESA_USE_CL
    cl::Buffer _offset_buf;
    #endif

    /**
     * @brief _kickdirection direction of the offset du to the kick
     */
    const Axis _kickdirection;

    /**
     * @brief _meshsize_kd size of the mesh in direction of the kick
     */
    #ifdef INOVESA_USE_CL
    const cl_int _meshsize_kd;
    #else
    const meshindex_t _meshsize_kd;
    #endif

    /**
     * @brief _meshsize_pd size of the mesh perpendicular to the kick
     */
    #ifdef INOVESA_USE_CL
    const cl_int _meshsize_pd;
    #else
    const meshindex_t _meshsize_kd;
    #endif

    /**
     * @brief updateHM
     *
     * does nothing when OpenCL is used
     */
    void updateHM();
};

}

#endif // KICKMAP_HPP
