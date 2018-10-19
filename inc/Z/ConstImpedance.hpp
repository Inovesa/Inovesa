// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
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
