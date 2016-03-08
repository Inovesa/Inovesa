#include "HM/RFKickMap.hpp"


vfps::RFKickMap::RFKickMap(vfps::PhaseSpace *in, vfps::PhaseSpace *out,
                           const vfps::meshindex_t xsize,
                           const vfps::meshindex_t ysize,
                           const vfps::meshaxis_t angle,
                           const vfps::HeritageMap::InterpolationType it)
    :
      KickMap(in,out,xsize,ysize,it)
{
    for(meshindex_t x=0; x<_xsize; x++) {
        _offset[x] = std::tan(angle)*(int(x)-int(_xsize/2));
    }
    updateHM();
}
