// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "SM/SourceMap.hpp"

#include <random>

namespace vfps
{

/**
 * @brief The FokkerPlanckMap class
 *
 * Source information for the FokkerPlanckMap is one dimensional
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

    ~FokkerPlanckMap() noexcept override;


    /**
     * @brief apply custom apply method needed to handle one dimensional SM
     */
    void apply() override;

    void applyTo(PhaseSpace::Position& pos) const override;

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
