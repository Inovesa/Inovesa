// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "Z/Impedance.hpp"

namespace vfps
{

class ParallelPlatesCSR : public Impedance
{
public:
    ParallelPlatesCSR( const size_t nfreqs
                     , const frequency_t f0
                     , const frequency_t f_max
                     , const double g
                     , oclhptr_t oclh = nullptr
                     );

private:
    static std::vector<vfps::impedance_t>
    __calcImpedance(const size_t nfreqs,
                    const frequency_t f0,
                    const frequency_t f_max,
                    const double g);
};

} // namespace VFPS
