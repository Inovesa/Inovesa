// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#include "Z/CollimatorImpedance.hpp"

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

vfps::CollimatorImpedance::CollimatorImpedance( const size_t n
                                              , const frequency_t f_max
                                              , const double outer
                                              , const double inner
                                              , oclhptr_t oclh
                                              )
    : ConstImpedance( n
                    , f_max
                    , { static_cast<frequency_t>(Z0/pi<double>()*std::log(outer/inner)),
                        static_cast<frequency_t>(0)}
                    , oclh
                    )
{
}
