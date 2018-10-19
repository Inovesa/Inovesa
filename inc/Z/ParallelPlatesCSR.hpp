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

class ParallelPlatesCSR : public Impedance
{
public:
    ParallelPlatesCSR( const size_t n
                     , const frequency_t f0
                     , const frequency_t f_max
                     , const double g
                     , oclhptr_t oclh = nullptr
                     );

private:
    static std::vector<vfps::impedance_t>
    __calcImpedance(const size_t n,
                    const frequency_t f0,
                    const frequency_t f_max,
                    const double g);
};

} // namespace VFPS
