#ifndef HERITAGEMAP_HPP
#define HERITAGEMAP_HPP

#include "PhaseSpace.hpp"

namespace vfps
{

class HeritageMap
{
public:
	/**
	 * @brief The INTERPOL_TYPE enum
	 */
	enum INTERPOL_TYPE
	{
		NONE=1,
		LINEAR=2,
		QUADRATIC=3,
		CUBIC=4
	};
	static constexpr INTERPOL_TYPE it = INTERPOL_TYPE::CUBIC;

public:
	HeritageMap(meshdata_t* in, meshdata_t* out, size_t xsize, size_t ysize);

	~HeritageMap();

	void apply();

protected:
	meshdata_t* const _data1D;
	meshdata_t* const _data1D_rotated;

	std::array<hi,it*it>** _heritage_map;
	std::array<hi,it*it>* const _heritage_map1D;

	const size_t _size;
	const size_t _xsize;
	const size_t _ysize;
};

}

#endif // HERITAGEMAP_HPP
