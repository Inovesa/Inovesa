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

#include "Z/ImpedanceFactory.hpp"

#include "IO/Display.hpp"

#include "Z/Impedance.hpp"
#include "Z/ConstImpedance.hpp"
#include "Z/CollimatorImpedance.hpp"
#include "Z/FreeSpaceCSR.hpp"
#include "Z/ParallelPlatesCSR.hpp"
#include "Z/ResistiveWall.hpp"

#include <boost/math/constants/constants.hpp>
using boost::math::constants::two_pi;

std::unique_ptr<vfps::Impedance>
vfps::makeImpedance( const size_t nfreqs
                   , oclhptr_t oclh
                   , const frequency_t fmax
                   , const double R_bend
                   , const double frev
                   , const double gap
                   , const bool use_csr
                   , const double s
                   , const double xi
                   , const double inner_coll_radius
                   , const std::string impedance_file
                   )
{
    auto f0 = physcons::c/(two_pi<double>()*R_bend);

    /*
     * Will create a default (zero) impedance and add different contributions.
     * If the impedance is still zero in the end (tracked by impedance_changed),
     * it will be replaced by a nullprt which is returned instead.
     */
    auto rv = std::make_unique<Impedance>( nfreqs,fmax, oclh);
    auto impedance_changed = false;

    if (gap != 0) {
        if (use_csr) {
            impedance_changed = true;
            Display::printText("... using CSR impedance");
            if (gap>0) {
                Display::printText("... shielded by parallel plates.");
                *rv += ParallelPlatesCSR(nfreqs,f0,fmax,gap);
            } else {
                Display::printText("... in free space.");
                *rv += FreeSpaceCSR(nfreqs,f0,fmax);
            }
        }
        auto radius = std::abs(gap/2);
        if ( s > 0 && xi >= -1 ) {
            impedance_changed = true;
            Display::printText("... using resistive wall impedance.");
            *rv += ResistiveWall(nfreqs,frev,fmax,physcons::c/frev,s,xi,radius);
        }
        if (0 < inner_coll_radius && inner_coll_radius < radius) {
            impedance_changed = true;
            Display::printText("... using collimator impedance.");
            *rv += CollimatorImpedance(nfreqs,fmax,radius,inner_coll_radius);
        }
    }

    if (impedance_file != "") {
        impedance_changed = true;
        Display::printText("Reading impedance from: \""
                           +impedance_file+"\"");
        *rv += Impedance(impedance_file,fmax);
    }

    // if impedance is still zero, a nullprt will be returned instead
    if (!impedance_changed) {
        Display::printText("... no impedance is used.");
        rv = nullptr;
    }

    return rv;
}

