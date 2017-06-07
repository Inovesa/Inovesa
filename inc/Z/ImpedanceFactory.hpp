#ifndef IMPEDANCEFACTORY_HPP
#define IMPEDANCEFACTORY_HPP

#include <memory>

#include "Z/Impedance.hpp"

namespace vfps
{

std::unique_ptr<Impedance> makeImpedance(const size_t nfreqs,
                                         const frequency_t fmax,
                                         const double frev,
                                         const double gap,
                                         const bool use_csr = true,
                                         const double s = 0,
                                         const double xi = 0,
                                         const double inner_coll_radius = 0,
                                         const std::string impedance_file = ""
                                         );

} // namespace vfps

#endif // IMPEDANCEFACTORY_HPP
