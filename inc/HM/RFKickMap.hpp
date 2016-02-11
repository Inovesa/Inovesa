#ifndef RFKICKMAP_HPP
#define RFKICKMAP_HPP

#include "HM/KickMap.hpp"

namespace vfps
{

class RFKickMap : public KickMap
{
public:
    RFKickMap(PhaseSpace* in, PhaseSpace* out,
              const meshindex_t xsize, const meshindex_t ysize,
              const meshaxis_t angle, const InterpolationType it);
};

} // namespace vfps

#endif // RFKICKMAP_HPP
