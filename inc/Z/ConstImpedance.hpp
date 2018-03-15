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

#ifndef CONSTIMPEDANCE_HPP
#define CONSTIMPEDANCE_HPP

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
                  #ifdef INOVESA_USE_OPENCL
                  , std::shared_ptr<OCLH> oclh = nullptr
                  #endif // INOVESA_USE_OPENCL
                  );

private:
    static std::vector<vfps::impedance_t>
    __calcImpedance(const size_t n,
                    const vfps::impedance_t Z);
};

} // namespace vfps

#endif // CONSTIMPEDANCE_HPP
