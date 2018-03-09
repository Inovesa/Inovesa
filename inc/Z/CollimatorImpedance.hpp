/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2017: Patrik Schönfeldt                                      *
 * Copyright (c) 2017: Karlsruhe Institute of Technology                      *
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

#ifndef COLLIMATORIMPEDANCE_HPP
#define COLLIMATORIMPEDANCE_HPP

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
                       , std::shared_ptr<OCLH> oclh = nullptr
                       );
};

} // namespace vfps

#endif // COLLIMATORIMPEDANCE_HPP
