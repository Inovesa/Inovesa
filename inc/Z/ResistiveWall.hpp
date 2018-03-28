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

#ifndef RESISTIVEWALL_HPP
#define RESISTIVEWALL_HPP

#include "Z/Impedance.hpp"

namespace vfps
{

/**
 * @brief The ResistiveWall class models the resistive wall impedance
 */
class ResistiveWall : public Impedance
{
public:
    /**
     * @brief ResistiveWall
     * @param n
     * @param f0
     * @param f_max
     * @param L: Length of the beam pipe
     * @param s: conductivity [S/m]
     * @param xi
     * @param b
     */
    ResistiveWall( const size_t n
                 , const frequency_t f0
                 , const frequency_t f_max
                 , const double L
                 , const double s
                 , const double xi
                 , const double b
                 , oclhptr_t oclh = nullptr
                 );

private:
    static std::vector<vfps::impedance_t>
    __calcImpedance(const size_t n,
                    const frequency_t f0,
                    const frequency_t f_max,
                    const double L,
                    const double s,
                    const double xi,
                    const double b);
};

} // namespace vfps

#endif // RESISTIVEWALL_HPP
