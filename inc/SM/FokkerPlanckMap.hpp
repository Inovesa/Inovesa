/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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

#ifndef FOKKERPLANCKMAP_HPP
#define FOKKERPLANCKMAP_HPP

#include "SM/SourceMap.hpp"

#include <random>

namespace vfps
{

/**
 * @brief The FokkerPlanckMap class
 *
 * Heritage information for the FokkerPlanckMap is one dimensional
 *
 * @todo Padding for faster memory access
 */
class FokkerPlanckMap : public SourceMap
{
public:
    /**
     * @brief The FPType enum holds ways to solve Fokker Planck equation
     */
    enum class FPType : uint_fast8_t {
        none=0,
        damping_only=1,
        diffusion_only=2,
        full=3
    };

    enum class FPTracking : uint_fast8_t {
        none=0,
        approximation1=1,
        approximation2=2,
        stochastic=3
    };

    enum DerivationType : uint_fast8_t {
        two_sided = 3,    // based on quadratic interpolation
        cubic = 4        // based on cubic interpolation
    };

public:
    FokkerPlanckMap() = delete;

    /**
     * @brief FokkerPlanckMap the (only) constructor
     * @param fpt any combination of damping/diffusion
     * @param e1 Marit: (deltat*2./(omegas*td))
     */
    FokkerPlanckMap(std::shared_ptr<PhaseSpace> in
                   , std::shared_ptr<PhaseSpace> out
                   , const meshindex_t xsize
                   , const meshindex_t ysize
                   , FPType fptype
                   , FPTracking fptrack, timeaxis_t e1
                   , DerivationType dt
                   , oclhptr_t oclh
                   );

    ~FokkerPlanckMap() noexcept;


    /**
     * @brief apply custom apply method needed to handle one dimensional SM
     */
    void apply() override;

    PhaseSpace::Position apply(PhaseSpace::Position pos) const override;

private:
    /**
     * @brief _dampincr damping decrement
     */
    timeaxis_t _dampdecr;

    mutable std::mt19937 _prng;

    mutable std::normal_distribution<meshaxis_t> _normdist ;

    const FPTracking _fptrack;

    const FPType _fptype;

    /**
     * @brief _meshxsize horizontal size of the meshaxis_t
     *
     * As this HeritageMap itself is one dimensional,
     * it has to save the size of the actual mesh seperatly.
     */
    const meshindex_t _meshxsize;
};

}

#endif // FOKKERPLANCKMAP_HPP
