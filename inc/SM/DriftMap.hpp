// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "SM/KickMap.hpp"

namespace vfps
{

class DriftMap : public KickMap
{
public:
    DriftMap(std::shared_ptr<PhaseSpace> in
            , std::shared_ptr<PhaseSpace> out
            , const meshindex_t xsize
            , const meshindex_t ysize
            , const std::vector<meshaxis_t> slip
            , const meshaxis_t E0
            , const InterpolationType it
            , const bool interpol_clamp
            , oclhptr_t oclh
            );

    #if INOVESA_ENABLE_CLPROFILING == 1
    ~DriftMap() noexcept override;
    #else
    ~DriftMap() = default;
    #endif // INOVESA_ENABLE_CLPROFILING
};

} // namespace fvps
