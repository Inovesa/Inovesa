#ifndef IMPEDANCEFACTORY_HPP
#define IMPEDANCEFACTORY_HPP

#include <memory>

#include "Z/Impedance.hpp"

namespace vfps
{

/**
 * @brief makeImpedance is a factory function for all kinds of impedances
 * @param nfreqs
 * @param fmax
 * @param f0
 * @param frev
 * @param gap
 * @param use_csr
 * @param s
 * @param xi
 * @param inner_coll_radius
 * @param impedance_file
 * @return pointer to fully initialized Impedance
 *
 * In priciple, an Impedance is easy to define and handle.
 * This factory function is designed with mainainability in mind:
 * It should be the single point in the program where Impedances are created.
 * As such, it provides a shortcut to common cases.
 */
std::unique_ptr<Impedance> makeImpedance( const size_t nfreqs
                                        #ifdef INOVESA_USE_OPENCL
                                        , std::shared_ptr<OCLH> oclh
                                        #endif // INOVESA_USE_OPENCL
                                        , const frequency_t fmax
                                        , const double f0
                                        , const double frev
                                        , const double gap
                                        , const bool use_csr = true
                                        , const double s = 0
                                        , const double xi = 0
                                        , const double inner_coll_radius = 0
                                        , const std::string impedance_file = ""
                                        );

} // namespace vfps

#endif // IMPEDANCEFACTORY_HPP
