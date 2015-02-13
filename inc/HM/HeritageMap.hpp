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
	HeritageMap(PhaseSpace* in, PhaseSpace* out, size_t xsize, size_t ysize);

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

#ifdef FR_USE_CL
	cl::Buffer _heritage_map1D_buf;
	cl::Kernel applyHM;

	cl::Buffer* data_read_buf;
	cl::Buffer* data_write_buf;

protected:
	void __initOpenCL();
#endif
};

}

#endif // HERITAGEMAP_HPP
