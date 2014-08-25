#ifndef MESH2D_HPP
#define MESH2D_HPP

#include <algorithm>
#include <array>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "CL/CLProgs.hpp"
#include "CL/OpenCLHandler.hpp"
#include "Ruler.hpp"

typedef float meshdata_t;

namespace vfps
{

typedef struct {
	unsigned int index;
	float weight;
} hi;

template <class data_t>
class Mesh2D
{	
public:
	/* choose number of meshpoints for interpolation;
	 * other values than 4 might not work properly because
	 * hardcoded dependencies are not yet flexible
	 */
	static constexpr unsigned int is = 4;

public:
	Mesh2D(std::array<Ruler<data_t>,2> axis) :
		_axis(axis),
		_data1D(new data_t[size<0>()*size<1>()]()),
		_data1D_rotated(new data_t[size<0>()*size<1>()]()),
		_heritage_map1D(new std::array<hi,is*is>[size<0>()*size<1>()]())
	{
		_data = new data_t*[size<0>()];
		_data_rotated = new data_t*[size<0>()];
		_heritage_map = new std::array<hi,16>*[size<0>()];
		for (unsigned int i=0; i<size<0>(); i++) {
			_data[i] = &(_data1D[i*size<1>()]);
			_data_rotated[i] = &(_data1D_rotated[i*size<1>()]);
			_heritage_map[i] = &(_heritage_map1D[i*size<1>()]);
		}
		_projection[0] = new data_t[size<0>()];
		_projection[1] = new data_t[size<1>()];
		_ws[0] = new data_t[size<0>()];
		_ws[1] = new data_t[size<1>()];

		const data_t ca = 3.;
		double dc = 1.;

		const data_t h03 = getDelta<0>()/3.;
		_ws[0][0] = h03;
		for (unsigned int x=1; x< size<0>()-1; x++){
			_ws[0][x] = h03 * (ca+dc);
			dc = -dc;
		}
		_ws[0][size<0>()-1] = h03;


		const data_t h13 = getDelta<1>()/3.;
		_ws[1][0] = h13;
		for (unsigned int x=1; x< size<1>()-1; x++){
			_ws[1][x] = h13 * (ca+dc);
			dc = -dc;
		}
		_ws[1][size<0>()-1] = h13;
    }

	Mesh2D(Ruler<data_t> axis1, Ruler<data_t> axis2) :
		Mesh2D(std::array<Ruler<data_t>,2>{{axis1,axis2}})
	{}

	Mesh2D(unsigned int xdim, data_t xmin, data_t xmax,
		   unsigned int ydim, data_t ymin, data_t ymax) :
		Mesh2D(Ruler<data_t>(xdim,xmin,xmax), Ruler<data_t>(ydim,ymin,ymax))
	{}

	Mesh2D(const Mesh2D& other) :
		Mesh2D(other._axis)
	{ std::copy_n(other._data1D,size<0>()*size<1>(),_data1D); }

	~Mesh2D()
	{
		delete [] _data;
		delete [] _data_rotated;
		delete [] _data1D;
		delete [] _data1D_rotated;
		delete [] _heritage_map;
		delete [] _heritage_map1D;
		delete [] _projection[0];
		delete [] _projection[1];
		delete [] _ws[0];
		delete [] _ws[1];
	}

    /**
     * @brief getData gives direct access to held data
     *
     * @return pointer to array holding size<0>()*size<1>() data points
     */
	inline const data_t* getData() const
	{
		data_synced.wait();
		return _data1D;
	}

	inline void syncData()
	{
		OCLH::queue.enqueueReadBuffer
				(_data1D_buf, CL_FALSE, 0,sizeof(float)*size<0>()*size<1>(),
				 _data1D,nullptr,&data_synced);
	}

    template <unsigned int x>
    inline const data_t getDelta() const
	{ return _axis[x].getDelta(); }

    template <unsigned int x>
    inline const data_t getMax() const
	{ return _axis[x].getMax(); }

    template <unsigned int x>
    inline const data_t getMin() const
	{ return _axis[x].getMin(); }

	template<unsigned int axis>
	data_t average()
	{
		if (axis == 0) {
			projectionToX();
		} else {
			projectionToY();
		}

		if (_moment[axis].size() == 0) {
			_moment[axis].resize(1);
		}

		data_t avg = 0;
		for (unsigned int i=0; i<size<axis>(); i++) {
			avg += i*_projection[axis][i];
		}
		avg = avg/size<axis>();

		_moment[axis][0] = avg;

		return x<axis>(avg);
	}

	template<unsigned int axis>
	data_t variance()
	{
		if (axis == 0) {
			projectionToX();
		} else {
			projectionToY();
		}

		if (_moment[axis].size() < 2) {
			_moment[axis].resize(2);
		}

		average();

		data_t avg = _moment[axis][0];
		data_t var = 0;
		for (unsigned int i=0; i<size<axis>(); i++) {
			var += (i-avg)*_projection[axis][i];
		}
		var = var/size<axis>();

		_moment[axis][1] = var;

		return x<axis>(var);
	}

	data_t* projectionToX() {
		for (unsigned int x=0; x < size<0>(); x++) {
			_projection[0][x] = 0;
			for (unsigned int y=0; y< size<1>(); y++) {
				_projection[0][x] += _data[x][y]*_ws[0][y];
			}
		}
		return _projection[0];
	}

	data_t* projectionToY() {
		for (unsigned int y=0; y< size<1>(); y++) {
			_projection[1][y] = 0;
			for (unsigned int x=0; x < size<0>(); x++) {
				_projection[1][y] += _data[x][y]*_ws[1][x];
			}
		}
		return _projection[1];
	}

    inline data_t* operator[](const unsigned int i) const
	{ return _data[i]; }

	Mesh2D& operator=(const Mesh2D& other)
    {
        if (_axis[0] == other._axis[0] && _axis[1] == other._axis[1]) {
            std::copy(other._data1D,other._data1D+size<0>()*size<1>(),_data1D);
        } else {
            throw std::range_error("Tried to assign two Mesh2D objects"
                                   "with different dimensions.");
        }
        return *this;
	}

	void setRotationMap(const data_t deltat)
	{
		constexpr double o3 = 1./3.;

		const double cos_dt = cos(deltat);
		const double sin_dt = sin(deltat);

		for(unsigned int p_i=0; p_i< size<1>(); p_i++) {
			for (unsigned int q_i=0; q_i< size<0>(); q_i++) {
				//  Find cell of inverse image (qp,pp) of grid point i,j.
				data_t qp = cos_dt*x<0>(q_i) - sin_dt*x<1>(p_i); //q', backward mapping
				data_t pp = sin_dt*x<0>(q_i) + cos_dt*x<1>(p_i); //p'
				//if ( willStayInMesh(qp,pp) )
				if ( insideMesh(qp,pp) )
				{
					/* choose number of meshpoints for interpolation;
					 * other values than 4 might not work properly because
					 * hardcoded dependencies are not yet flexible
					 */
					constexpr unsigned int interpolation_steps = 4;

					unsigned int id = floor((qp-getMin<0>())/getDelta<0>()); //meshpoint smaller q'
					unsigned int jd = floor((pp-getMin<1>())/getDelta<1>()); //numper of lower mesh point from p'

					// arrays of Lagrange interpolation
					std::array<data_t,interpolation_steps> laq;
					std::array<data_t,interpolation_steps> lap;

					// gridpoint matrix used for interpolation
					std::array<std::array<hi,interpolation_steps>,interpolation_steps> ph;

					//Scaled arguments of interpolation functions:
					data_t xiq = (qp-x<0>(id))/getDelta<0>();  //distance from id
					data_t xip = (pp-x<1>(jd))/getDelta<1>(); //distance of p' from lower mesh point

					for (unsigned int j1=0; j1<interpolation_steps; j1++) {
						unsigned int j0 = jd+j1-1;
						for (unsigned int i1=0; i1<interpolation_steps; i1++) {
							unsigned int i0 = id+i1-1;
							if(i0< size<0>() && j0 < size<1>() ){
								ph[i1][j1].index = i0*size<1>()+j0;
							} else {
								ph[i1][j1].index = 0;
							}
						}
					}

					//  Vectors of Lagrange interpolation functions, not including factors of 1/2.
					laq[0] = (xiq-1)*(xiq-2);
					laq[1] = (xiq+1)*laq[0];
					laq[0] *= -xiq*o3;
					laq[3] = xiq*(xiq+1);
					laq[2] = -(xiq-2)*laq[3];
					laq[3] *= (xiq-1)*o3;

					lap[0] = (xip-1)*(xip-2);
					lap[1] = (xip+1)*lap[0];
					lap[0] *= -xip*o3;
					lap[3] = xip*(xip+1);
					lap[2] = -(xip-2)*lap[3];
					lap[3] *= (xip-1)*o3;

					//  Assemble Lagrange interpolation as quadratic form, restoring factors of 1/2:
					for (size_t i1=0; i1<interpolation_steps; i1++) {
						for (size_t j1=0; j1<interpolation_steps; j1++){
							(ph[i1][j1]).weight = laq[i1] * lap[j1] * 0.25;
							_heritage_map[q_i][p_i][i1*interpolation_steps+j1]
									= ph[i1][j1];
						}
					}
				}
			}
		}
	}

	/**
	 * @brief rotate
	 * @param deltat
	 */
	void rotate()
	{

		#define FR_USE_CL
		#ifdef FR_USE_CL
		OCLH::queue.enqueueNDRangeKernel (
				applyHM1D,
				cl::NullRange,
				cl::NDRange(size<0>()*size<1>()));
		#ifdef CL_VERSION_1_2
		OCLH::queue.enqueueBarrierWithWaitList();
		#else // CL_VERSION_1_2
		OCLH::queue.enqueueBarrier();
		#endif // CL_VERSION_1_2
		OCLH::queue.enqueueCopyBuffer(_data1D_rotated_buf,
									  _data1D_buf,
									  0,0,
									  sizeof(float)*size<0>()*size<1>());
		#else // FR_USE_CL
		for (unsigned int i=0; i< size<0>()*size<1>(); i++) {
			_data1D_rotated[i] = 0.0;
			for (hi h: _heritage_map1D[i]) {
				_data1D_rotated[i] += _data1D[h.index]*h.weight;
			}
		}
		std::copy_n(_data1D_rotated,size<0>()*size<1>(),_data1D);
		#endif // FR_USE_CL
	}

	void kick(const std::vector<data_t> &AF)
	{
		;
	}

	/**
	 * @brief rotateAndKick
	 * @param deltat
	 * @param AF
	 *
	 * Method has to be reviewed!
	 */
	void rotateAndKick(const data_t deltat, const std::vector<data_t> &AF)
	{
		constexpr double o3 = 1./3.;

		const double cos_dt = cos(deltat);
		const double sin_dt = sin(deltat);

		for(unsigned int p_i=0; p_i< size<1>(); p_i++) {
			for (unsigned int q_i=0; q_i< size<0>(); q_i++) {
				//  Find cell of inverse image (qp,pp) of grid point i,j.
				data_t qp = cos_dt*x<0>(q_i) - sin_dt*(x<1>(p_i)+AF[q_i]); //q', backward mapping
				data_t pp = sin_dt*x<0>(q_i) + cos_dt*(x<1>(p_i)+AF[q_i]); //p'
				//if (willStayInMesh(qp,pp)) {
				if (insideMesh(qp,pp)) {
					unsigned int id = floor((qp-getMin<0>())/getDelta<0>()); //meshpoint smaller q'
					unsigned int jd = floor((pp-getMin<1>())/getDelta<1>()); //numper of lower mesh point from p'

					// arrays of Lagrange interpolation
					std::array<data_t,is> laq;
					std::array<data_t,is> lap;

					// gridpoint matrix used for interpolation
					std::array<std::array<data_t,is>,is> ph;

					//Scaled arguments of interpolation functions:
					data_t xiq = (qp-x<0>(id))/getDelta<0>();  //distance from id
					data_t xip = (pp-x<1>(jd))/getDelta<1>(); //distance of p' from lower mesh point

					for (unsigned int j1=0; j1<is; j1++) {
						unsigned int j0 = jd+j1-1;
						for (unsigned int i1=0; i1<is; i1++) {
							unsigned int i0 = id+i1-1;
							if(i0< size<0>() && j0 < size<1>() ){
								ph[i1][j1] = _data[i0][j0];
							} else {
								ph[i1][j1]=0.0;
							}
						}
					}

					//  Vectors of Lagrange interpolation functions, not including factors of 1/2.
					laq[0] = (xiq-1)*(xiq-2);
					laq[1] = (xiq+1)*laq[0];
					laq[0] *= -xiq*o3;
					laq[3] = xiq*(xiq+1);
					laq[2] = -(xiq-2)*laq[3];
					laq[3] *= (xiq-1)*o3;

					lap[0] = (xip-1)*(xip-2);
					lap[1] = (xip+1)*lap[0];
					lap[0] *= -xip*o3;
					lap[3] = xip*(xip+1);
					lap[2] = -(xip-2)*lap[3];
					lap[3] *= (xip-1)*o3;

					//  Assemble Lagrange interpolation as quadratic form, restoring factors of 1/2:
					_data_rotated[q_i][p_i] = 0.;
					for (size_t i1=0; i1<is; i1++) {
						for (size_t j1=0; j1<is; j1++){
							_data_rotated[q_i][p_i] += ph[i1][j1] * laq[i1] * lap[j1];
						}
					}
					_data_rotated[q_i][p_i] *= 0.25;
				} else {
					_data_rotated[q_i][p_i] = 0.;
				}
			}
		}
		swapDataTmp();
	}

	void swapDataTmp() {
		data_t** data_swap = _data;
		_data = _data_rotated;
		_data_rotated = data_swap;
	}

    template <unsigned int x>
    inline const unsigned int size() const
	{ return _axis[x].getNSteps(); }

	template <unsigned int axis>
	inline const data_t x(const unsigned int n) const
	{ return _axis[axis][n]; }

protected:
	const std::array<Ruler<data_t>,2> _axis;

	std::array<data_t*,2> _projection;

	data_t** _data;

	data_t** _data_rotated;

	data_t* const _data1D;

	data_t* const _data1D_rotated;

	cl::Event data_synced;

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
	std::array<std::vector<data_t>,2> _moment;

	std::array<data_t*,2> _ws;

	inline bool insideMesh(data_t x, data_t y) const
	{
		return (x <= getMax<0>() && y <= getMax<1>() &&
				x >= getMin<0>() && y >= getMin<1>() );
	}

	inline bool willStayInMesh(data_t x, data_t y) const
	{
		return (sqrt(pow(x,2)+pow(y,2)) < getMax<0>());
	}

private:
	cl::Buffer _data1D_buf;
	cl::Buffer _data1D_rotated_buf;
	cl::Buffer _heritage_map1D_buf;
	cl::Kernel applyHM1D;

public:
	void __initOpenCL()
	{
		_data1D_buf = cl::Buffer(OCLH::context,
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				sizeof(float)*size<0>()*size<1>(), _data1D);
		_data1D_rotated_buf = cl::Buffer(OCLH::context,
				CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR,
				sizeof(float)*size<0>()*size<1>(), _data1D_rotated);
		_heritage_map1D_buf = cl::Buffer(OCLH::context,
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				sizeof(hi)*is*is*size<0>()*size<1>(),
				_heritage_map1D);
		applyHM1D = cl::Kernel(CLProgApplyHM::p, "applyHM1D");
		applyHM1D.setArg(0, _data1D_buf);
		applyHM1D.setArg(1, _heritage_map1D_buf);
		applyHM1D.setArg(2, _data1D_rotated_buf);
	}
};

}

#endif // MESH2D_HPP
