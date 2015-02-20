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
#include "defines.hpp"
#include "Ruler.hpp"

namespace vfps
{

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
	inline meshdata_t* getData() const
	{
		return _data1D;
	}

	inline meshaxis_t getDelta(const unsigned int x) const
	{ return _axis[x].getDelta(); }

	inline meshaxis_t getMax(const unsigned int x) const
	{ return _axis[x].getMax(); }

	inline meshaxis_t getMin(const unsigned int x) const
	{ return _axis[x].getMin(); }

	meshdata_t average(const unsigned int axis);

	meshdata_t variance(const unsigned int axis);

	integral_t* projectionToX();

	integral_t* projectionToY();

	inline meshdata_t* operator[](const unsigned int i) const
	{ return _data[i]; }

	PhaseSpace& operator=(PhaseSpace other);

	inline unsigned int nMeshCells() const
	{ return _axis[0].getNSteps()*_axis[1].getNSteps(); }

	inline unsigned int nMeshCells(const unsigned int x) const
	{ return _axis[x].getNSteps(); }

	inline meshaxis_t size(const unsigned int x) const
	{ return _axis[x].size(); }

	inline meshaxis_t x(const unsigned int axis,
							  const unsigned int n) const
	{ return _axis[axis][n]; }

	/**
	 * @brief swap
	 * @param other
	 */
	friend void swap(PhaseSpace& first, PhaseSpace& second) noexcept;

protected:
	const std::array<Ruler<meshaxis_t>,2> _axis;

	std::array<integral_t*,2> _projection;

	meshdata_t** _data;

	meshdata_t* _data1D;

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

#ifdef INOVESA_USE_CL
public:
	cl::Buffer data_buf;
#endif
};

/**
 * @brief swap
 * @param first
 * @param second
 *
 * @todo adjust to also swap cl::Buffer
 */
void swap(PhaseSpace& first, PhaseSpace& second) noexcept;

}

#endif // PHASESPACE_HPP
