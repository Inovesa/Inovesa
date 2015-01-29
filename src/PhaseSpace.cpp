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
	double dc = 1.;

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
		avg += i*_projection[axis][i];
	}
	avg = avg/size(axis);

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
		var += (i-avg)*_projection[axis][i];
	}
	var = var/size(axis);

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

	meshaxis_t cos_dt = cos(deltat);
	meshaxis_t sin_dt = sin(deltat);

	std::string interpol_str;
	switch (mn) {
	case ROTATION_TYPE::MESH:
		interpol_str = "mesh";
		break;
	case ROTATION_TYPE::NORMAL:
		interpol_str = "normal";
		break;
	case ROTATION_TYPE::SPACE:
	default:
		interpol_str = "space";
		break;
	}

	std::ofstream hm("hermasum_"+interpol_str+".dat");
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
				xiq = std::modf(qp, &qq_int);
				xip = std::modf(pp, &qp_int);
				id = qq_int;
				jd = qp_int;
				break;
			case ROTATION_TYPE::NORMAL:
				qp = cos_dt*((q_i-size(0)/2.0)/size(0))
						- sin_dt*((p_i-size(1)/2.0)/size(1))+0.5;
				pp = sin_dt*((q_i-size(0)/2.0)/size(0))
						+ cos_dt*((p_i-size(1)/2.0)/size(1))+0.5;
				qcoord = qp*size(0);
				pcoord = pp*size(1);
				xiq = std::modf(qcoord, &qq_int);
				xip = std::modf(pcoord, &qp_int);
				id = qq_int;
				jd = qp_int;
				break;
			case ROTATION_TYPE::SPACE:
				qp = cos_dt*x(0,q_i) - sin_dt*x(1,p_i) - getMin(0);
				pp = sin_dt*x(0,q_i) + cos_dt*x(1,p_i) - getMin(1);
				qcoord = qp/getDelta(0);
				pcoord = pp/getDelta(1);
				xiq = std::modf(qcoord, &qq_int);
				xip = std::modf(pcoord, &qp_int);
				id = qq_int;
				jd = qp_int;
				break;
			}

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
					icq[0] = (1-xiq);
					icq[1] = xiq;

					icp[0] = (1-xip);
					icp[1] = xip;
					break;

				case INTERPOL_TYPE::QUADRATIC:
					icq[0] = xiq*(xiq-1)/2;
					icq[1] = 1-xiq*xiq;
					icq[2] = xiq*(xiq+1)/2;

					icp[0] = xip*(xip-1)/2;
					icp[1] = 1-xip*xip;
					icp[2] = xip*(xip+1)/2;

					break;

				case INTERPOL_TYPE::CUBIC:
					constexpr double o3 = 1./3.;

					icq[0] = (xiq-1)*(xiq-2)/2;
					icq[1] = (xiq+1)*icq[0];
					icq[0] *= -o3*xiq;
					icq[3] = xiq*(xiq+1)/2;
					icq[2] = (2-xiq)*icq[3];
					icq[3] *= (xiq-1)*o3;

					icp[0] = (xip-1)*(xip-2)/2;
					icp[1] = (xip+1)*icp[0];
					icp[0] *= -o3*xip;
					icp[3] = xip*(xip+1)/2;
					icp[2] = (2-xip)*icp[3];
					icp[3] *= (xip-1)*o3;
					break;
				}

				//  Assemble interpolation
				for (size_t hmq=0; hmq<it; hmq++) {
					for (size_t hmp=0; hmp<it; hmp++){
						hmc[hmp*it+hmq] = icq[hmp]*icp[hmq];
					}
				}

				// renormlize to minimize rounding errors
				//normalize(it,hmc.data());

				// write heritage map
				for (unsigned int j1=0; j1<it; j1++) {
					unsigned int j0 = jd+j1-1;
					for (unsigned int i1=0; i1<it; i1++) {
						unsigned int i0 = id+i1-1;
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
			interpol_t hmsum = 0;
			for (hi h : _heritage_map[q_i][p_i]) {
				hmsum += h.weight;
				ti[h.index] += h.weight;
			}

			long int epsilon;
			if (std::is_same<vfps::interpol_t,double>::value) {
				epsilon = std::round((static_cast<double>(hmsum)
									  -1.0+DBL_EPSILON)/DBL_EPSILON)-1.0;
			}
			if (std::is_same<vfps::interpol_t,float>::value) {
				epsilon = std::round((static_cast<float>(hmsum)
									  -1.0f+FLT_EPSILON)/FLT_EPSILON)-1.0f;
			}
			if (std::is_same<vfps::interpol_t,Share>::value) {
				epsilon = static_cast<long int>(hmsum);
			}

			if (p_i > 1 && q_i > 1
				&& p_i < size(1)-1
				&& q_i < size(0) -1 ) {
				if (epsilon >= -5 && epsilon <= 5) {
					num[epsilon+5]++;
				} else {
					num[11]++;
				}
			}
			hm << epsilon << '\t';
		}
		hm << std::endl;
	}
	std::ofstream hmb("hermabin_"+interpol_str+".dat");
	for (unsigned int i=0; i< num.size(); i++) {
		hmb << int(i)-5 << '\t' << num[i] << std::endl;
	}
	std::ofstream tm("target_"+interpol_str+".dat");
	for (unsigned int x=0; x<size(0); x++) {
		for (unsigned int y=0; y<size(1); y++) {
			tm << static_cast<float>(ti[y*size(0)+x] - 1.0) << '\t';
		}
		tm << std::endl;
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
			_data1D_rotated[i] += _data1D[h.index]*h.weight;
		}
		// handle overflow
		if (std::is_same<vfps::meshdata_t,unsigned int>::value) {
			if (_data1D_rotated[i] > UINT32_MAX/4*3) {
				_data1D_rotated[i] = 0;
			}
		}
	}
	std::copy_n(_data1D_rotated,size(0)*size(1),_data1D);
#endif // FR_USE_CL
}

void vfps::PhaseSpace::kick(const std::vector<meshdata_t>& AF)
{
	;
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
