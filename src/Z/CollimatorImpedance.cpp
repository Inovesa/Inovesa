/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2017-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2017-2018: Karlsruhe Institute of Technology                 *
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

#include "Z/CollimatorImpedance.hpp"

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

vfps::CollimatorImpedance::CollimatorImpedance(const size_t n,
                                   const frequency_t f_max,
                                   const double outer,
                                   const double inner)
    :
      ConstImpedance( n,f_max,
          {static_cast<frequency_t>(Z0/pi<double>()*std::log(outer/inner)),
           static_cast<frequency_t>(0)}
      )
{
}
