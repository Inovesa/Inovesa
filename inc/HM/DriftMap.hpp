#ifndef DRIFTMAP_HPP
#define DRIFTMAP_HPP

#include "HM/KickMap.hpp"

namespace vfps
{

class DriftMap : public KickMap
{
public:
	DriftMap(PhaseSpace* in, PhaseSpace* out,
			 const meshindex_t xsize, const meshindex_t ysize,
			 const meshaxis_t angle, const InterpolationType it);
};

} // namespace fvps

#endif // DRIFTMAP_HPP
