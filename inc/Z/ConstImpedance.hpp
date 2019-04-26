// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "Z/Impedance.hpp"

namespace vfps
{

/**
 * @brief The CollimatorImpedance class models
 * the impedance effect of a collimator
 */
class ConstImpedance : public Impedance
{
public:
    ConstImpedance( const size_t n
                  , const frequency_t f_max
                  , const vfps::impedance_t Z
                  , oclhptr_t oclh = nullptr
                  );

private:
    static std::vector<vfps::impedance_t>
    __calcImpedance(const size_t n,
                    const vfps::impedance_t Z);
};

} // namespace vfps
