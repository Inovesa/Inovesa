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

#ifndef PHASESPACEFACTORY_HPP
#define PHASESPACEFACTORY_HPP

#include <memory>
#include <string>

#include "PS/PhaseSpace.hpp"

namespace vfps {

#ifdef INOVESA_USE_HDF5
std::unique_ptr<PhaseSpace> makePSFromHDF5(std::string fname,
                                           int64_t startdiststep,
                                           meshaxis_t qmin, meshaxis_t qmax,
                                           meshaxis_t pmin, meshaxis_t pmax,
                                           const double bunch_charge,
                                           const double bunch_current,
                                           double xscale, double yscale);
#endif // INOVESA_USE_HDF5

std::unique_ptr<PhaseSpace> makePSFromPNG(std::string fname,
                                          meshaxis_t qmin, meshaxis_t qmax,
                                          meshaxis_t pmin, meshaxis_t pmax,
                                          const double bunch_charge,
                                          const double bunch_current,
                                          double xscale, double yscale);

std::unique_ptr<PhaseSpace> makePSFromTXT(std::string fname, int64_t ps_size,
                                          meshaxis_t qmin, meshaxis_t qmax,
                                          meshaxis_t pmin, meshaxis_t pmax,
                                          const double bunch_charge,
                                          const double bunch_current,
                                          double xscale, double yscale);


} // namespace vfps

#endif // PHASESPACEFACTORY_HPP
