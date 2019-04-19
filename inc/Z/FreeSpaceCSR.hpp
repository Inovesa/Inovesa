// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "Z/Impedance.hpp"

namespace vfps
{

class FreeSpaceCSR : public Impedance
{
public:
    FreeSpaceCSR( const size_t n
                , const frequency_t f_rev
                , const frequency_t f_max
                , oclhptr_t oclh = nullptr
                );

private:
    static std::vector<vfps::impedance_t>
    __calcImpedance(const size_t n,
                    const frequency_t f_rev,
                    const frequency_t f_max);
};

} // namespace vfps
