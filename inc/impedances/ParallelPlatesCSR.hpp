#ifndef PARALLELPLATESCSR_HPP
#define PARALLELPLATESCSR_HPP

#include "impedances/Impedance.hpp"

namespace vfps
{

class ParallelPlatesCSR : public Impedance
{
public:
    ParallelPlatesCSR(const size_t n,
                      const frequency_t f_rev,
                      const frequency_t f_max,
                      const double h);

private:
    static std::vector<vfps::impedance_t>
    __calcImpedance(const size_t n,
                    const frequency_t f_rev,
                    const frequency_t f_max,
                    const double h);
};

} // namespace VFPS

#endif // PARALLELPLATESCSR_HPP
