#include "PhaseSpace.hpp"

vfps::PhaseSpace::PhaseSpace(std::array<Ruler<meshaxis_t>,2> axis) :
	_axis(axis),
	_data1D(new meshdata_t[size(0)*size(1)]()),
	_data1D_rotated(new meshdata_t[size(0)*size(1)]()),
	_heritage_map1D(new std::array<hi,it*it>[size(0)*size(1)]())
{
	_data = new meshdata_t*[size(0)];
	_data_rotated = new meshdata_t*[size(0)];
	_heritage_map = new std::array<hi,it*it>*[size(0)];
	for (unsigned int i=0; i<size(0); i++) {
		_data[i] = &(_data1D[i*size(1)]);
		_data_rotated[i] = &(_data1D_rotated[i*size(1)]);
		_heritage_map[i] = &(_heritage_map1D[i*size(1)]);
	}
	_projection[0] = new meshdata_t[size(0)];
	_projection[1] = new meshdata_t[size(1)];
	_ws[0] = new meshdata_t[size(0)];
	_ws[1] = new meshdata_t[size(1)];

	const meshdata_t ca = 3.;
	meshdata_t dc = 1;

	const meshdata_t h03 = getDelta(0)/3.;
	_ws[0][0] = h03;
	for (unsigned int x=1; x< size(0)-1; x++){
		_ws[0][x] = h03 * (ca+dc);
		dc = -dc;
	}
	_ws[0][size(0)-1] = h03;


	const meshdata_t h13 = getDelta(1)/3.;
	_ws[1][0] = h13;
	for (unsigned int x=1; x< size(1)-1; x++){
		_ws[1][x] = h13 * (ca+dc);
		dc = -dc;
	}
	_ws[1][size(0)-1] = h13;

	img_size = {{size(0),size(1)}};
}


vfps::PhaseSpace::PhaseSpace(Ruler<meshaxis_t> axis1, Ruler<meshaxis_t> axis2) :
	PhaseSpace(std::array<Ruler<meshaxis_t>,2>{{axis1,axis2}})
{}

vfps::PhaseSpace::PhaseSpace(unsigned int xdim, meshaxis_t xmin, meshaxis_t xmax,
							 unsigned int ydim, meshaxis_t ymin, meshaxis_t ymax) :
	PhaseSpace(Ruler<meshaxis_t>(xdim,xmin,xmax),
			   Ruler<meshaxis_t>(ydim,ymin,ymax))
{}

vfps::PhaseSpace::PhaseSpace(const vfps::PhaseSpace& other) :
	PhaseSpace(other._axis)
{ std::copy_n(other._data1D,size(0)*size(1),_data1D); }

vfps::PhaseSpace::~PhaseSpace()
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

vfps::meshdata_t vfps::PhaseSpace::average(const unsigned int axis)
{
	if (axis == 0) {
		projectionToX();
	} else {
		projectionToY();
	}

	if (_moment[axis].size() == 0) {
		_moment[axis].resize(1);
	}

	meshdata_t avg = 0;
	for (unsigned int i=0; i<size(axis); i++) {
		avg += meshdata_t(i)*_projection[axis][i];
	}
	avg = avg/meshdata_t(size(axis));

	_moment[axis][0] = avg;

	return x(axis,avg);
}

vfps::meshdata_t vfps::PhaseSpace::variance(const unsigned int axis)
{
	if (axis == 0) {
		projectionToX();
	} else {
		projectionToY();
	}

	if (_moment[axis].size() < 2) {
		_moment[axis].resize(2);
	}

	average(axis);

	meshdata_t avg = _moment[axis][0];
	meshdata_t var = 0;
	for (unsigned int i=0; i<size(axis); i++) {
		var += (meshdata_t(i)-avg)*_projection[axis][i];
	}
	var = var/meshdata_t(size(axis));

	_moment[axis][1] = var;

	return x(axis,var);
}

vfps::meshdata_t* vfps::PhaseSpace::projectionToX() {
	for (unsigned int x=0; x < size(0); x++) {
		_projection[0][x] = 0;
		for (unsigned int y=0; y< size(1); y++) {
			_projection[0][x] += _data[x][y]*_ws[0][y];
		}
	}
	return _projection[0];
}

vfps::meshdata_t* vfps::PhaseSpace::projectionToY() {
	for (unsigned int y=0; y< size(1); y++) {
		_projection[1][y] = 0;
		for (unsigned int x=0; x < size(0); x++) {
			_projection[1][y] += _data[x][y]*_ws[1][x];
		}
	}
	return _projection[1];
}

vfps::PhaseSpace& vfps::PhaseSpace::operator=(const vfps::PhaseSpace& other)
{
	if (_axis[0] == other._axis[0] && _axis[1] == other._axis[1]) {
		std::copy(other._data1D,other._data1D+size(0)*size(1),_data1D);
	} else {
		throw std::range_error("Tried to assign two Mesh2D objects"
							   "with different dimensions.");
	}
	return *this;
}

void vfps::PhaseSpace::setRotationMap(const meshaxis_t deltat,
									  const vfps::PhaseSpace::ROTATION_TYPE mn)
{
	std::vector<interpol_t> ti;
	ti.resize(size(0)*size(1));

	const meshaxis_t cos_dt = cos(deltat);
	const meshaxis_t sin_dt = sin(deltat);

	std::array<unsigned int,12> num;
	num.fill(0);
	for(unsigned int p_i=0; p_i< size(1); p_i++) {
		for (unsigned int q_i=0; q_i< size(0); q_i++) {

			// Cell of inverse image (qp,pp) of grid point i,j.
			meshaxis_t qp; //q', backward mapping
			meshaxis_t pp; //p'
			// interpolation type specific q and p coordinates
			meshaxis_t pcoord;
			meshaxis_t qcoord;
			meshaxis_t qq_int;
			meshaxis_t qp_int;
			//Scaled arguments of interpolation functions:
			unsigned int id; //meshpoint smaller q'
			unsigned int jd; //numper of lower mesh point from p'
			interpol_t xiq; //distance from id
			interpol_t xip; //distance of p' from lower mesh point
			switch (mn) {
			case ROTATION_TYPE::MESH:
				qp = cos_dt*(q_i-size(0)/2.0)
						- sin_dt*(p_i-size(1)/2.0)+size(0)/2.0;
				pp = sin_dt*(q_i-size(0)/2.0)
						+ cos_dt*(p_i-size(1)/2.0)+size(1)/2.0;
				qcoord = qp;
				pcoord = pp;
				break;
			case ROTATION_TYPE::NORMAL:
				qp = cos_dt*((q_i-size(0)/2.0)/size(0))
						- sin_dt*((p_i-size(1)/2.0)/size(1));
				pp = sin_dt*((q_i-size(0)/2.0)/size(0))
						+ cos_dt*((p_i-size(1)/2.0)/size(1));
				qcoord = (qp+0.5)*size(0);
				pcoord = (pp+0.5)*size(1);
				break;
			case ROTATION_TYPE::NORMAL2:
				qp = cos_dt*(2*static_cast<int>(q_i)-static_cast<int>(size(0)))/static_cast<int>(size(0))
				   - sin_dt*(2*static_cast<int>(p_i)-static_cast<int>(size(1)))/static_cast<int>(size(1));

				pp = sin_dt*(2*static_cast<int>(q_i)-static_cast<int>(size(0)))/static_cast<int>(size(0))
				   + cos_dt*(2*static_cast<int>(p_i)-static_cast<int>(size(1)))/static_cast<int>(size(1));
				qcoord = (qp+1)*size(0)/2;
				pcoord = (pp+1)*size(1)/2;
				break;
			case ROTATION_TYPE::SPACE:
				qp = cos_dt*x(0,q_i) - sin_dt*x(1,p_i);
				pp = sin_dt*x(0,q_i) + cos_dt*x(1,p_i);
				qcoord = (qp-getMin(0))/getDelta(0);
				pcoord = (pp-getMin(1))/getDelta(1);
				break;
			}
			xiq = std::modf(qcoord, &qq_int);
			xip = std::modf(pcoord, &qp_int);
			id = qq_int;
			jd = qp_int;

			if (id <  size(0) && jd < size(1))
			{
				// gridpoint matrix used for interpolation
				std::array<std::array<hi,it>,it> ph;

				// arrays of interpolation coefficients
				std::array<interpol_t,it> icq;
				std::array<interpol_t,it> icp;

				std::array<interpol_t,it*it> hmc;

				// create vectors containing interpolation coefficiants
				switch (it)
				{
				case INTERPOL_TYPE::NONE:
					icq[0] = 1;

					icp[0] = 1;
					break;

				case INTERPOL_TYPE::LINEAR:
					icq[0] = interpol_t(1)-xiq;
					icq[1] = xiq;

					icp[0] = interpol_t(1)-xip;
					icp[1] = xip;
					break;

				case INTERPOL_TYPE::QUADRATIC:
					icq[0] = xiq*(xiq-interpol_t(1))/interpol_t(2);
					icq[1] = interpol_t(1)-xiq*xiq;
					icq[2] = xiq*(xiq+interpol_t(1))/interpol_t(2);

					icp[0] = xip*(xip-interpol_t(1))/interpol_t(2);
					icp[1] = interpol_t(1)-xip*xip;
					icp[2] = xip*(xip+interpol_t(1))/interpol_t(2);
					break;

				case INTERPOL_TYPE::CUBIC:
					icq[0] = (xiq-interpol_t(1))*(xiq-interpol_t(2))*xiq
							* interpol_t(-1./6.);
					icq[1] = (xiq+interpol_t(1))*(xiq-interpol_t(1))
							* (xiq-interpol_t(2)) / interpol_t(2);
					icq[2] = (interpol_t(2)-xiq)*xiq*(xiq+interpol_t(1))
							/ interpol_t(2);
					icq[3] = xiq*(xiq+interpol_t(1))*(xiq-interpol_t(1))
							* interpol_t(1./6.);

					icp[0] = (xip-interpol_t(1))*(xip-interpol_t(2))*xip
							* interpol_t(-1./6.);
					icp[1] = (xip+interpol_t(1))*(xip-interpol_t(1))
							* (xip-interpol_t(2)) / interpol_t(2);
					icp[2] = (interpol_t(2)-xip)*xip*(xip+interpol_t(1))
							/ interpol_t(2);
					icp[3] = xip*(xip+interpol_t(1))*(xip-interpol_t(1))
							* interpol_t(1./6.);
					break;
				}

				//  Assemble interpolation
				for (size_t hmq=0; hmq<it; hmq++) {
					for (size_t hmp=0; hmp<it; hmp++){
						hmc[hmp*it+hmq] = icq[hmp]*icp[hmq];
					}
				}


				// renormlize to minimize rounding errors
//				renormalize(hmc.size(),hmc.data());

				// write heritage map
				for (unsigned int j1=0; j1<it; j1++) {
					unsigned int j0 = jd+j1-(it-1)/2;
					for (unsigned int i1=0; i1<it; i1++) {
						unsigned int i0 = id+i1-(it-1)/2;
						if(i0< size(0) && j0 < size(1) ){
							ph[i1][j1].index = i0*size(1)+j0;
							ph[i1][j1].index2d = {{i0,j0}};
							ph[i1][j1].weight = hmc[i1*it+j1];
						} else {
							ph[i1][j1].index = 0;
							ph[i1][j1].index2d = {{0,0}};
							ph[i1][j1].weight = 0;
						}
						_heritage_map[q_i][p_i][i1*it+j1]
								= ph[i1][j1];
					}
				}
			}
		}
	}
}

void vfps::PhaseSpace::rotate()
{
#ifdef FR_USE_CL
	OCLH::queue.enqueueNDRangeKernel (
				applyHM2D,
				cl::NullRange,
				cl::NDRange(size(0),size(1)));
#ifdef CL_VERSION_1_2
	OCLH::queue.enqueueBarrierWithWaitList();
#else // CL_VERSION_1_2
	OCLH::queue.enqueueBarrier();
#endif // CL_VERSION_1_2
	cl::size_t<3> null3d;
	cl::size_t<3> imgsize;
	imgsize[0] = size(0);
	imgsize[1] = size(1);
	imgsize[2] = 1;
	OCLH::queue.enqueueCopyImage(_data_rotated_buf,
								 _data_buf,
								 null3d,null3d,imgsize);
#else // FR_USE_CL
	for (unsigned int i=0; i< size(0)*size(1); i++) {
		_data1D_rotated[i] = 0;
		for (hi h: _heritage_map1D[i]) {
			_data1D_rotated[i] += _data1D[h.index]*static_cast<meshdata_t>(h.weight);
		}
		// handle overshooting
		meshdata_t ceil=std::numeric_limits<fixp32>::min();
		meshdata_t flor=std::numeric_limits<fixp32>::max();
		for (size_t x=1; x<=2; x++) {
			for (size_t y=1; y<=2; y++) {
				ceil = std::max(ceil,_data1D[_heritage_map1D[i][x*it+y].index]);
				flor = std::min(flor,_data1D[_heritage_map1D[i][x*it+y].index]);
			}
		}
		_data1D_rotated[i] = std::max(std::min(ceil,_data1D_rotated[i]),flor);
	}
	std::copy_n(_data1D_rotated,size(0)*size(1),_data1D);
#endif // FR_USE_CL
}

void vfps::PhaseSpace::kick(const std::vector<meshdata_t>& AF)
{
	;
}

void vfps::PhaseSpace::renormalize(size_t n, float* args)
{
	float sum=0;
	for( size_t i=0; i<n; i++) {
		sum += args[i];
	}

	if (sum > 0) {
		for( size_t i=0; i<n; i++) {
			args[i] /= sum;
		}
	} else {
		for( size_t i=0; i<n; i++) {
			args[i] = 1.0f/n;
		}
	}
}

void vfps::PhaseSpace::renormalize(size_t n, fixp32* args)
{
	std::vector<fixp32> tmpshare;
	tmpshare.reserve(n);
	fixp32 sum=0;
	for( size_t i=0; i<n; i++) {
		sum += args[i];
		tmpshare.push_back(args[i]);
	}

	if (sum == fixp32(0)) {
		for( size_t i=0; i<n; i++) {
			tmpshare[i] = 1;
		}
		sum = n;
	}
	std::list<std::pair<fixp32,size_t>> remainders;
	fixp32 rensum=0;
	for( size_t i=0; i<n; i++) {
		fixp32 val = tmpshare[i];
		fixp32 remainder = val - sum*(val/sum);
		val /= sum;
		rensum += val;
		args[i] = val;
		remainders.push_back(std::pair<fixp32,size_t>(remainder,i));
	}

	fixp32 rest = fixp32(1) - rensum;
	if (rest > fixp32(0)) {
		// sort remainders descending
		remainders.sort([](	const std::pair<fixp32,size_t> &lhs,
							const std::pair<fixp32,size_t> &rhs)
						-> bool
						{ return lhs.first > rhs.first; }
		);
		do {
			args[remainders.front().second] += std::numeric_limits<fixp32>::epsilon();
			remainders.pop_front();
			rest -= fixp32(std::numeric_limits<fixp32>::epsilon());
		} while (rest > fixp32(0));
	} else  if (rest < fixp32(0)) {
		// sort remainders ascending
		remainders.sort([](	const std::pair<fixp32,size_t> &lhs,
							const std::pair<fixp32,size_t> &rhs)
						-> bool
						{ return lhs.first < rhs.first; }
		);
		do {
			args[remainders.front().second] -= std::numeric_limits<fixp32>::epsilon();
			remainders.pop_front();
			rest += fixp32(std::numeric_limits<fixp32>::epsilon());
		} while (rest > fixp32(0));
	}
}

void vfps::PhaseSpace::__initOpenCL()
{
	_data_buf = cl::Image2D(OCLH::context,
							CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
							cl::ImageFormat(CL_R,CL_FLOAT),
							size(0),size(1),0, _data1D);
	_data_rotated_buf = cl::Image2D(OCLH::context,
									CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR,
									cl::ImageFormat(CL_R,CL_FLOAT),
									size(0),size(1),0, _data1D_rotated);
	_heritage_map1D_buf = cl::Buffer(OCLH::context,
									 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
									 sizeof(hi)*it*it*size(0)*size(1),
									 _heritage_map1D);
	applyHM2D = cl::Kernel(CLProgApplyHM::p, "applyHM2D");
	applyHM2D.setArg(0, _data_buf);
	applyHM2D.setArg(1, _heritage_map1D_buf);
	applyHM2D.setArg(2, size(1));
	applyHM2D.setArg(3, it*it);
	applyHM2D.setArg(4, _data_rotated_buf);
}
