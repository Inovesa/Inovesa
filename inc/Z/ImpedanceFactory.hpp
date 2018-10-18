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

#pragma once

#include <memory>

#include "Z/Impedance.hpp"

namespace vfps
{

/**
 * @brief makeImpedance is a factory function for all kinds of impedances
 * @param nfreqs
 * @param fmax
 * @param f0
 * @param frev
 * @param gap
 * @param use_csr
 * @param s
 * @param xi
 * @param inner_coll_radius
 * @param impedance_file
 * @return pointer to fully initialized Impedance
 *
 * In priciple, an Impedance is easy to define and handle.
 * This factory function is designed with mainainability in mind:
 * It should be the single point in the program where Impedances are created.
 * As such, it provides a shortcut to common cases.
 */
std::unique_ptr<Impedance> makeImpedance(const size_t nfreqs
                                        , oclhptr_t oclh
                                        , const frequency_t fmax
                                        , const double R_bend
                                        , const double frev
                                        , const double gap
                                        , const bool use_csr = true
                                        , const double s = 0
                                        , const double xi = 0
                                        , const double inner_coll_radius = 0
                                        , const std::string impedance_file = ""
                                        );

} // namespace vfps
