#ifndef PHASESPACE_HPP
#define PHASESPACE_HPP

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <GL/gl.h>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "CL/CLProgs.hpp"
#include "CL/OpenCLHandler.hpp"
#include "Ruler.hpp"

#define FR_USE_GUI
#define FR_USE_CL
#define FR_CL_SYNC_BLOCKING CL_TRUE

typedef float meshdata_t;
typedef float interpol_t;

namespace vfps
{

typedef struct {
	unsigned int index;
	cl_uint2 index2d;
	interpol_t weight;
} hi;

class PhaseSpace
{	
public:
	/**
	 * @brief The MESH_NORMALIZATION enum
	 */
	enum class ROTATION_TYPE
	{
		MESH=0, NORMAL, SPACE
	};

	enum INTERPOL_TYPE
	{
		NONE=1,
		LINEAR=2,
		QUADRATIC=3,
		CUBIC=4
	};
	static constexpr INTERPOL_TYPE is = INTERPOL_TYPE::CUBIC;

public:
	PhaseSpace(std::array<Ruler<interpol_t>,2> axis);

	PhaseSpace(Ruler<interpol_t> axis1, Ruler<interpol_t> axis2);

	PhaseSpace(	unsigned int xdim, interpol_t xmin, interpol_t xmax,
				unsigned int ydim, interpol_t ymin, interpol_t ymax);

	PhaseSpace(const PhaseSpace& other);

	~PhaseSpace();

	 /**
	  * @brief getData gives direct access to held data
	  *
	  * @return pointer to array holding size<0>()*size<1>() data points
	  */
	inline const meshdata_t* getData() const
	{
#ifdef FR_USE_CL
#if (FR_CL_SYNC_BLOCKING == CL_FALSE)
		data_synced.wait();
#endif // FR_CL_SYNC_BLOCKING
#endif // FR_USE_CL
		return _data1D;
	}

	inline void syncData()
	{
		cl::size_t<3> null3d;
		cl::size_t<3> imgsize;
		imgsize[0] = size(0);
		imgsize[1] = size(1);
		imgsize[2] = 1;
		OCLH::queue.enqueueReadImage
					(_data_buf, FR_CL_SYNC_BLOCKING, null3d,imgsize,0,0,
					_data1D,nullptr,&data_synced);
	}

	inline const interpol_t getDelta(const unsigned int x) const
	{ return _axis[x].getDelta(); }

	inline const interpol_t getMax(const unsigned int x) const
	{ return _axis[x].getMax(); }

	inline const interpol_t getMin(const unsigned int x) const
	{ return _axis[x].getMin(); }

	meshdata_t average(const unsigned int axis);

	meshdata_t variance(const unsigned int axis);

	meshdata_t* projectionToX();

	meshdata_t* projectionToY();

	inline meshdata_t* operator[](const unsigned int i) const
	{ return _data[i]; }

	PhaseSpace& operator=(const PhaseSpace& other);

	void setRotationMap(const interpol_t deltat,
						const ROTATION_TYPE mn=ROTATION_TYPE::SPACE);

	/**
	 * @brief rotate applies the rotation map
	 */
	void rotate();

	void kick(const std::vector<meshdata_t> &AF);

	inline const unsigned int size(const unsigned int x) const
	{ return _axis[x].getNSteps(); }

	inline const interpol_t x(const unsigned int axis,
							  const unsigned int n) const
	{ return _axis[axis][n]; }

protected:
	const std::array<Ruler<interpol_t>,2> _axis;

	std::array<meshdata_t*,2> _projection;

	meshdata_t** _data;

	meshdata_t** _data_rotated;

	meshdata_t* const _data1D;

	meshdata_t* const _data1D_rotated;

	cl::Event data_synced;

	cl_uint2 img_size;

	std::array<hi,is*is>** _heritage_map;
	std::array<hi,is*is>* const _heritage_map1D;

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

	inline bool insideMesh(interpol_t x, interpol_t y) const
	{
		return (x <= getMax(0) && y <= getMax(1) &&
				x >= getMin(0) && y >= getMin(1) );
	}

	inline bool willStayInMesh(interpol_t x, interpol_t y) const
	{
		return (sqrt(pow(x,2)+pow(y,2)) < getMax(0));
	}

private:
	cl::Image2D _data_buf;
	cl::Image2D _data_rotated_buf;
	cl::Buffer _heritage_map1D_buf;
	cl::Kernel applyHM1D;
	cl::Kernel applyHM2D;
	cl::Kernel rotateImg;

public:
	void __initOpenCL();
};

}

#endif // PHASESPACE_HPP
