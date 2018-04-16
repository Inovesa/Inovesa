/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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

#include "MessageStrings.hpp"

const std::string vfps::copyright_notice() noexcept {
    std::string rv (
        "Inovesa Numerical Optimized Vlasov-Equation Solver Application\n"
        "\n");
    if (
            #if FXP_FRACPART < 31
            std::is_same<vfps::meshdata_t,fixp32>::value ||
            #endif
            std::is_same<vfps::meshdata_t,fixp64>::value)
    {
        rv += "Copyright (c) 2007-2009 Peter Schregle (FPML)\n";
    }
    rv+="Copyright (c) 2012-2018 Patrik Schönfeldt\n"
        "Copyright (c) 2014-2018 Karlsruhe Institute of Technology\n"
        "Copyright (c) 2017-2018 Johannes Schestag\n"
        "Copyright (c) 2017 Tobias Boltz\n"
        "Copyright (c) 2017 Patrick Schreiber\n"
        "Copyright (c) 2017 Julian Gethmann\n"
        "Copyright (c) 2017 Matthias Blaicher\n"
        "Copyright (c) 1997-2016 John C. Bowman,\n"
        "\tUniversity of Alberta (Array class)\n"
        #ifdef INOVESA_USE_OPENGL
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
        if ( INOVESA_VERSION_FIX >= 0 ) {
            sstream << '.' << INOVESA_VERSION_FIX;
        } else if ( INOVESA_VERSION_FIX == -1 ) {
                sstream << " alpha";
        } else if ( INOVESA_VERSION_FIX == -2 ) {
            sstream << " beta";
        } else {
            sstream << " RC" << std::abs(INOVESA_VERSION_FIX)-2;
        }
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
                #ifdef DEBUG
                << std::endl << '\t' << "DEBUG"
                #endif
                #ifdef INOVESA_ENABLE_INTERRUPT
                << std::endl << '\t' << "INOVESA_ENABLE_INTERRUPT"
                #endif
                #ifdef INOVESA_ENABLE_CLPROFILING
                << std::endl << '\t' << "INOVESA_ENABLE_CLPROFILING"
                #endif
                #ifdef INOVESA_USE_OPENCL
                << std::endl << '\t' << "INOVESA_USE_OPENCL"
                #endif
                #ifdef INOVESA_USE_CLFFT
                << std::endl << '\t' << "INOVESA_USE_CLFFT"
                #endif
                #ifdef INOVESA_USE_OPENGL
                << std::endl << '\t' << "INOVESA_USE_OPENGL"
                #endif
                #ifdef INOVESA_USE_HDF5
                << std::endl << '\t' << "INOVESA_USE_HDF5"
                #endif
                #ifdef INOVESA_USE_PNG
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

