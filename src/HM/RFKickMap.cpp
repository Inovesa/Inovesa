#include "inc/HM/RFKickMap.hpp"


vfps::RFKickMap::RFKickMap(vfps::PhaseSpace *in, vfps::PhaseSpace *out,
                           const vfps::meshindex_t xsize,
                           const vfps::meshindex_t ysize,
                           const vfps::meshaxis_t angle,
                           const vfps::HeritageMap::InterpolationType it)
    :
      KickMap(in,out,xsize,ysize,it)
{
    ;
}
