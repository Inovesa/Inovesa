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

#ifndef FREESPACECSR_HPP
#define FREESPACECSR_HPP

#include "Z/Impedance.hpp"

namespace vfps
{

/**
 * @brief The ResistiveWall class models the resistive wall impedance
 *
 * First line according to Eq. 13.6 in
 * "Coupling Impedances and Beam Instabilities in Accelerator Rings"
 * by K.Y. Ng
 *
 * Z(w) = [1+j*sign(w)] sqrt(abs(w)*mu/2/s)*R/b
 *      = [1+j*sign(w)] sqrt(pi*mu*f_rev/s)*R/b*sqrt(abs(n))
 *
 * w: frequency ( = 2*pi*n*f_rev)
 * mu: magnetic permeability
 * s: conductivity
 * R: mean radius of the accelerator ring
 * b: inside radius of the beam pipe
 */
class ResistiveWall : public Impedance
{
public:
    ResistiveWall(const size_t n,
                  const frequency_t f_rev,
                  const frequency_t f_max,
                  const double mu,
                  const double s,
                  const double b);

private:
    static std::vector<vfps::impedance_t>
    __calcImpedance(const size_t n,
                    const frequency_t f_rev,
                    const frequency_t f_max,
                    const double mu,
                    const double s,
                    const double b);
};

} // namespace vfps

#endif // FREESPACECSR_HPP
