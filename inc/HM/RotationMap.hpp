#ifndef ROTATIONMAP_HPP
#define ROTATIONMAP_HPP

#include "HeritageMap.hpp"

namespace vfps
{

class RotationMap : public HeritageMap
{
public:
	RotationMap(PhaseSpace* in, PhaseSpace* out,
				const size_t xsize, const size_t ysize,
				const meshaxis_t angle);
};

}

#endif // ROTATIONMAP_HPP
