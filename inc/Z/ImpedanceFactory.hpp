#ifndef IMPEDANCEFACTORY_HPP
#define IMPEDANCEFACTORY_HPP

#include <memory>

#include "Z/Impedance.hpp"

namespace vfps
{

std::unique_ptr<Impedance> makeImpedance(const size_t nfreqs,
                                         const frequency_t fmax,
                                         const std::string impedance_file,
                                         const double frev,
                                         const bool use_csr,
                                         const double gap,
                                         const double s,
                                         const double xi,
                                         const double inner_coll_radius);

} // namespace vfps

#endif // IMPEDANCEFACTORY_HPP
