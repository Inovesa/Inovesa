// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include "Z/ConstImpedance.hpp"

namespace vfps
{

/**
 * @brief The CollimatorImpedance class models
 * the impedance effect of a collimator
 */
class CollimatorImpedance : public ConstImpedance
{
public:
    CollimatorImpedance( const size_t n
                       , const frequency_t f_max
                       , const double outer
                       , const double inner
                       , oclhptr_t oclh = nullptr
                       );
};

} // namespace vfps
