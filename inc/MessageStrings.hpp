/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2017: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2017: Karlsruhe Institute of Technology                 *
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

#ifndef MESSAGESTRINGS_HPP
#define MESSAGESTRINGS_HPP

#include <iomanip>
#include <memory>
#include <sstream>


#include "PS/PhaseSpace.hpp"

namespace vfps
{

const std::string copyright_notice() noexcept;

const std::string inovesa_version();

const std::string status_string(std::shared_ptr<PhaseSpace> ps,
                                float roatation,
                                float total_rotations);

#ifdef INOVESA_ENABLE_CLPROFILING
const std::string printProfilingInfo(const cl::vector<cl::Event> &events);
#endif

} // namespace vfps

#endif // MESSAGESTRINGS_HPP
