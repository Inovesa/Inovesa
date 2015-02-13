#ifndef PHASESPACE_HPP
#define PHASESPACE_HPP

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <GL/gl.h>
#include <list>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "CL/CLProgs.hpp"
#include "CL/OpenCLHandler.hpp"
#include "fixed_point.h"
#include "Ruler.hpp"

#define FR_USE_GUI
#define FR_USE_CL
#define FR_CL_SYNC_BLOCKING CL_TRUE


namespace vfps
{
typedef fpml::fixed_point<int,2,29> fixp32;

constexpr unsigned int ps_xsize = 512;
constexpr unsigned int ps_ysize = 512;

typedef float meshaxis_t;
typedef float meshdata_t;
typedef float interpol_t;

typedef struct {
	unsigned int index;
	cl_uint2 index2d;
	interpol_t weight;
} hi;

class PhaseSpace
{	
public:
	PhaseSpace(std::array<Ruler<meshaxis_t>,2> axis);

	PhaseSpace(Ruler<meshaxis_t> axis1, Ruler<meshaxis_t> axis2);

	PhaseSpace(	meshaxis_t xmin, meshaxis_t xmax,
				meshaxis_t ymin, meshaxis_t ymax);

	PhaseSpace(const PhaseSpace& other);

	~PhaseSpace();

	 /**
	  * @brief getData gives direct access to held data
	  *
	  * @return pointer to array holding size<0>()*size<1>() data points
	  */
	inline meshdata_t* const getData() const
	{
#ifdef FR_USE_CL
#if (FR_CL_SYNC_BLOCKING == CL_FALSE)
		data_synced.wait();
#endif // FR_CL_SYNC_BLOCKING
#endif // FR_USE_CL
		return _data1D;
	}

#ifdef FR_USE_CL
	inline void syncData()
	{
		OCLH::queue.enqueueReadBuffer
					(data_buf, FR_CL_SYNC_BLOCKING,
					 0,sizeof(float)*size(0)*size(1),
					_data1D,nullptr,&data_synced);
	}
#endif

	inline const meshaxis_t getDelta(const unsigned int x) const
	{ return _axis[x].getDelta(); }

	inline const meshaxis_t getMax(const unsigned int x) const
	{ return _axis[x].getMax(); }

	inline const meshaxis_t getMin(const unsigned int x) const
	{ return _axis[x].getMin(); }

	meshdata_t average(const unsigned int axis);

	meshdata_t variance(const unsigned int axis);

	meshdata_t* projectionToX();

	meshdata_t* projectionToY();

	inline meshdata_t* operator[](const unsigned int i) const
	{ return _data[i]; }

	PhaseSpace& operator=(PhaseSpace other);

	inline const unsigned int size(const unsigned int x) const
	{ return _axis[x].getNSteps(); }

	inline const meshaxis_t x(const unsigned int axis,
							  const unsigned int n) const
	{ return _axis[axis][n]; }

	/**
	 * @brief swap
	 * @param other
	 */
	friend void swap(PhaseSpace& first, PhaseSpace& second) noexcept;

protected:
	const std::array<Ruler<meshaxis_t>,2> _axis;

	std::array<meshdata_t*,2> _projection;

	meshdata_t** _data;

	meshdata_t* _data1D;

	cl::Event data_synced;

	cl_uint2 img_size;

	/**
	 * @brief _moment: holds the moments for distributions
	 *			in both axis in mesh coordinates
	 *
	 * 0: average
	 * 1: variance
	 * 2: skewness
	 * 3: kurtosis
	 */
	std::array<std::vector<meshdata_t>,2> _moment;

	std::array<meshdata_t*,2> _ws;

#ifdef FR_USE_CL
public:
	cl::Buffer data_buf;

public:
	void __initOpenCL();
#endif
};

void swap(PhaseSpace& first, PhaseSpace& second) noexcept;

}

#endif // PHASESPACE_HPP
