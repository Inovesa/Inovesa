#include "Z/ImpedanceFactory.hpp"

#include "IO/Display.hpp"

#include "Z/Impedance.hpp"
#include "Z/CollimatorImpedance.hpp"
#include "Z/FreeSpaceCSR.hpp"
#include "Z/ParallelPlatesCSR.hpp"
#include "Z/ResistiveWall.hpp"

std::unique_ptr<vfps::Impedance> vfps::makeImpedance(const size_t nfreqs,
                                               const frequency_t fmax,
                                               const std::string impedance_file,
                                               const double frev,
                                               const bool use_csr,
                                               const double gap,
                                               const double s,
                                               const double xi,
                                               const double inner_coll_radius)
{
    auto rv = std::make_unique<Impedance>(nfreqs,fmax);

    if (gap != 0) {
        if (use_csr) {
            Display::printText("Using CSR impedance");
            if (gap>0) {
                Display::printText("... shielded by parallel plates.");
                *rv += ParallelPlatesCSR(nfreqs,frev,fmax,gap);
            } else {
                Display::printText("... in free space.");
                *rv += FreeSpaceCSR(nfreqs,frev,fmax);
            }
        }
        auto radius = std::abs(gap/2);
        if ( s > 0 && xi >= -1 ) {
            Display::printText("... using resistive wall impedance.");
            *rv += ResistiveWall(nfreqs,frev,fmax, s,xi,radius);
        }
        if (inner_coll_radius > radius) {
            Display::printText("... using collimator impedance.");
            *rv += CollimatorImpedance(nfreqs,fmax,radius,inner_coll_radius);
        }
    }

    if (impedance_file != "") {
        Display::printText("Reading impedance from: \""
                           +impedance_file+"\"");
        *rv += Impedance(impedance_file,fmax);
    }

    return rv;
}

