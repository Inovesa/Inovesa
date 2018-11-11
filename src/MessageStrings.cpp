// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#include <boost/config.hpp>

#include "MessageStrings.hpp"


const std::string vfps::copyright_notice() noexcept {
    std::string rv (
        "Inovesa Numerical Optimized Vlasov-Equation Solver Application\n"
        "\n");
    rv+="Copyright (c) 2012-2018 Patrik Schönfeldt\n"
        "Copyright (c) 2014-2018 Karlsruhe Institute of Technology\n"
        "Copyright (c) 2017-2018 Johannes Schestag\n"
        "Copyright (c) 2017 Tobias Boltz\n"
        "Copyright (c) 2017 Patrick Schreiber\n"
        "Copyright (c) 2017 Julian Gethmann\n"
        "Copyright (c) 2017 Matthias Blaicher\n"
        "Copyright (c) 1997-2016 John C. Bowman,\n"
        "\tUniversity of Alberta (Array class)\n"
        #if INOVESA_USE_OPENGL == 1
        "Copyright (c) 2014-2015 Nathaniel J. Smith, Stefan van der Walt\n"
        "\t (Magma color code)\n"
        #endif
        "\n"
        "Inovesa is free software: you can redistribute it and/or modify\n"
        "it under the terms of the GNU General Public License as published by\n"
        "the Free Software Foundation, either version 3 of the License, or\n"
        "(at your option) any later version.\n"
        "\n"
        "Inovesa is distributed in the hope that it will be useful,\n"
        "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
        "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
        "GNU General Public License for more details.\n"
        "\n"
        "You should have received a copy of the GNU General Public License"
        "along with Inovesa. If not, see <http://www.gnu.org/licenses/>.";
    return rv;
}

const std::string vfps::inovesa_version(const bool verbose) {
    std::stringstream sstream;
    auto branchstr = std::string(GIT_BRANCH);
    if (branchstr.empty() || branchstr == "master" || branchstr.front()=='v') {
        // version number release branches (and releases)
        sstream << 'v' << INOVESA_VERSION_MAJOR
                << '.' << INOVESA_VERSION_MINOR;
        #if ( INOVESA_VERSION_FIX >= 0)
            sstream << '.' << INOVESA_VERSION_FIX;
        #elif ( INOVESA_VERSION_FIX == -1 )
                sstream << " alpha";
        #elif ( INOVESA_VERSION_FIX == -2 )
            sstream << " beta";
        #else
            sstream << " RC" << std::abs(INOVESA_VERSION_FIX)-2;
        #endif
        if (branchstr.front()=='v') {
            // plus commit for pre-release branches
            sstream << ", Commit: "<< GIT_COMMIT;
        }
     } else {
        // no version number for development or feature branches
        sstream << "Branch: "<< GIT_BRANCH
                << ", Commit: "<< GIT_COMMIT;
     }
    if (verbose) {
        sstream << std::endl << "Build options:"
                << std::endl << "Compiler:"
                << std::endl << '\t' << BOOST_COMPILER
                << std::endl << "Build options:"
                #if DEBUG == 1
                << std::endl << '\t' << "DEBUG"
                #endif
                #if INOVESA_ENABLE_INTERRUPT == 1
                << std::endl << '\t' << "INOVESA_ENABLE_INTERRUPT"
                #endif
                #if INOVESA_ENABLE_CLPROFILING == 1
                << std::endl << '\t' << "INOVESA_ENABLE_CLPROFILING"
                #endif
                #if INOVESA_USE_OPENCL == 1
                << std::endl << '\t' << "INOVESA_USE_OPENCL"
                #endif
                #if INOVESA_USE_CLFFT == 1
                << std::endl << '\t' << "INOVESA_USE_CLFFT"
                #endif
                #if INOVESA_USE_OPENGL == 1
                << std::endl << '\t' << "INOVESA_USE_OPENGL"
                #endif
                #if INOVESA_USE_HDF5 == 1
                << std::endl << '\t' << "INOVESA_USE_HDF5"
                #endif
                #if INOVESA_USE_PNG == 1
                << std::endl << '\t' << "INOVESA_USE_PNG"
                #endif
                ;
    }

    return sstream.str();
}



const std::string vfps::status_string(std::shared_ptr<PhaseSpace> ps,
                                      float roatation,
                                      float total_rotations)
{
    std::stringstream status;
    status.precision(5);
    status << std::setw(6) << roatation << '/' << total_rotations;
    float delta = 1.0 - ps->getIntegral();
    status << "\t1-Q/Q_0=";
    // this is to have constant spacing in the output
    if (delta < 0) {
        status << "-";
    } else {
        status << "+";
    }
    status << std::scientific << std::setw(10)
           << std::abs(delta)
           << "\ts_p="
           << std::fixed << std::setprecision(4)
           << *ps->getEnergySpread();

    return status.str();
}

