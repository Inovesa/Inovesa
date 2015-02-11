#ifndef ROTATIONMAP_HPP
#define ROTATIONMAP_HPP

#include "HeritageMap.hpp"

namespace vfps
{

class RotationMap : public HeritageMap
{
public:
	/**
	 * @brief The ROTATION_TYPE enum
	 */
	enum class ROTATION_TYPE
	{
		MESH=0, NORMAL, NORMAL2
	};

public:
	RotationMap(meshdata_t* in, meshdata_t* out,
				const size_t xsize, const size_t ysize,
				const meshaxis_t angle,
				const ROTATION_TYPE mn=ROTATION_TYPE::NORMAL2);
};

}

#endif // ROTATIONMAP_HPP
