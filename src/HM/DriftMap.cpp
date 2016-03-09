#include "inc/HM/DriftMap.hpp"

vfps::DriftMap::DriftMap(PhaseSpace *in, PhaseSpace *out,
						 const meshindex_t xsize,
						 const meshindex_t ysize,
						 const meshaxis_t angle,
						 const InterpolationType it)
	:
	  KickMap(in,out,xsize,ysize,it,DirectionOfKick::x)
{
	for(meshindex_t y=0; y<_ysize; y++) {
		_offset[y] = std::tan(angle)*( static_cast<int>(y)
									  -static_cast<int>(_ysize/2));
	}
	updateHM();
}
