#include "HM/RFKickMap.hpp"


vfps::RFKickMap::RFKickMap(PhaseSpace *in, PhaseSpace *out,
                           const meshindex_t xsize,
                           const meshindex_t ysize,
                           const meshaxis_t angle,
                           const InterpolationType it)
    :
      KickMap(in,out,xsize,ysize,it,DirectionOfKick::y)
{
    for(meshindex_t x=0; x<_xsize; x++) {
        _offset[x] = std::tan(angle)*( static_cast<int>(_xsize/2)
                                      -static_cast<int>(x));
    }
    updateHM();
}
