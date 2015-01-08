#include "PhaseSpace.hpp"

vfps::PhaseSpace::PhaseSpace(std::array<Ruler<interpol_t>,2> axis) :
	_axis(axis),
	_data1D(new meshdata_t[size(0)*size(1)]()),
	_data1D_rotated(new meshdata_t[size(0)*size(1)]()),
	_heritage_map1D(new std::array<hi,is*is>[size(0)*size(1)]())
{
	_data = new meshdata_t*[size(0)];
	_data_rotated = new meshdata_t*[size(0)];
	_heritage_map = new std::array<hi,16>*[size(0)];
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


vfps::PhaseSpace::PhaseSpace(Ruler<interpol_t> axis1, Ruler<interpol_t> axis2) :
	PhaseSpace(std::array<Ruler<interpol_t>,2>{{axis1,axis2}})
{}

vfps::PhaseSpace::PhaseSpace(unsigned int xdim, interpol_t xmin, interpol_t xmax,
							 unsigned int ydim, interpol_t ymin, interpol_t ymax) :
	PhaseSpace(Ruler<interpol_t>(xdim,xmin,xmax),
			   Ruler<interpol_t>(ydim,ymin,ymax))
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

meshdata_t* vfps::PhaseSpace::projectionToX() {
	for (unsigned int x=0; x < size(0); x++) {
		_projection[0][x] = 0;
		for (unsigned int y=0; y< size(1); y++) {
			_projection[0][x] += _data[x][y]*_ws[0][y];
		}
	}
	return _projection[0];
}

meshdata_t* vfps::PhaseSpace::projectionToY() {
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

void vfps::PhaseSpace::setRotationMap(const interpol_t deltat,
									  const vfps::PhaseSpace::ROTATION_TYPE mn)
{
	std::vector<interpol_t> ti;
	ti.resize(size(0)*size(1));

	constexpr double o3 = 1./3.;

	interpol_t cos_dt = cos(deltat);
	interpol_t sin_dt = sin(deltat);

#ifdef FR_USE_CL
	rotation = {{cos_dt,sin_dt}};
#endif

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
	std::array<unsigned int,8> num;
	num.fill(0);
	for(unsigned int p_i=0; p_i< size(1); p_i++) {
		for (unsigned int q_i=0; q_i< size(0); q_i++) {
			bool valid_meshpoint;

			// Cell of inverse image (qp,pp) of grid point i,j.
			interpol_t qp;
			interpol_t pp;
			// interpolation type specific q and p coordinates
			interpol_t pcoord;
			interpol_t qcoord;
			interpol_t qq_int;
			interpol_t qp_int;
			//Scaled arguments of interpolation functions:
			unsigned int id; //meshpoint smaller q'
			unsigned int jd; //numper of lower mesh point from p'
			interpol_t xiq; //distance from id
			interpol_t xip; //distance of p' from lower mesh point
			switch (mn) {
			case ROTATION_TYPE::MESH:
				qp = cos_dt*(q_i-size(0)/2.0)
						- sin_dt*(p_i-size(1)/2.0)+size(0)/2.0; //q', backward mapping
				pp = sin_dt*(q_i-size(0)/2.0)
						+ cos_dt*(p_i-size(1)/2.0)+size(1)/2.0; //p'
				xiq = std::modf(qp, &qq_int);
				id = qq_int; //meshpoint smaller q'
				xip = std::modf(pp, &qp_int);
				jd = qp_int; //numper of lower mesh point from p'

				break;
			case ROTATION_TYPE::NORMAL:
				qp = cos_dt*((q_i-size(0)/2.0)/size(0))
						- sin_dt*((p_i-size(1)/2.0)/size(1))+0.5; //q', backward mapping
				pp = sin_dt*((q_i-size(0)/2.0)/size(0))
						+ cos_dt*((p_i-size(1)/2.0)/size(1))+0.5; //p'
				qcoord = qp*size(0);
				pcoord = pp*size(1);
				xiq = std::modf(qcoord, &qq_int);
				xip = std::modf(pcoord, &qp_int);
				id = qq_int;
				jd = qp_int;
				break;
			case ROTATION_TYPE::SPACE:
			default:
				qp = cos_dt*x(0,q_i) - sin_dt*x(1,p_i); //q', backward mapping
				pp = sin_dt*x(0,q_i) + cos_dt*x(1,p_i); //p'
				valid_meshpoint =  insideMesh(qp,pp);
				if (valid_meshpoint) {
					id = floor((qp-getMin(0))/getDelta(0));
					jd = floor((pp-getMin(1))/getDelta(1));
					xiq = (qp-x(0,id))/getDelta(0);
					xip = (pp-x(1,jd))/getDelta(1);
				}
				break;
			}

			switch (mn) {
			case ROTATION_TYPE::MESH:
				valid_meshpoint =  (id <  size(0) && jd < size(1));
				break;
			case ROTATION_TYPE::NORMAL:
				valid_meshpoint =  ( qp > 0.0 && qp < 1.0 && pp > 0.0 && pp < 1.0 );
				break;
			default:
				break;
			}
			if (valid_meshpoint)
			{
				/* choose number of meshpoints for interpolation;
					 * other values than 4 might not work properly because
					 * hardcoded dependencies are not yet flexible
					 */
				constexpr unsigned int interpolation_steps = 4;

				// arrays of Lagrange interpolation
				std::array<interpol_t,interpolation_steps> laq;
				std::array<interpol_t,interpolation_steps> lap;

				// gridpoint matrix used for interpolation
				std::array<std::array<hi,interpolation_steps>,interpolation_steps> ph;

				for (unsigned int j1=0; j1<interpolation_steps; j1++) {
					unsigned int j0 = jd+j1-1;
					for (unsigned int i1=0; i1<interpolation_steps; i1++) {
						unsigned int i0 = id+i1-1;
						if(i0< size(0) && j0 < size(1) ){
							ph[i1][j1].index = i0*size(1)+j0;
							ph[i1][j1].index2d = {{i0,j0}};
						} else {
							ph[i1][j1].index = 0;
							ph[i1][j1].index2d = {{0,0}};
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
						(ph[i1][j1]).weight = laq[i1] * lap[j1] /4;
						_heritage_map[q_i][p_i][i1*interpolation_steps+j1]
								= ph[i1][j1];
					}
				}
			}
			double hmsum = 0.0;
			for (hi h : _heritage_map[q_i][p_i]) {
				hmsum += h.weight;
				ti[h.index] += h.weight;
			}
			int epsilon = std::round((hmsum-1.0+DBL_EPSILON)/DBL_EPSILON)-1;
			if (p_i > 1 && q_i > 1
				&& p_i < size(1)-1
				&& q_i < size(0) -1 ) {
				if (epsilon >= -3 && epsilon <= 3) {
					num[epsilon+3]++;
				} else {
					num[7]++;
				}
			}
			hm << epsilon << '\t';
		}
		hm << std::endl;
	}
	std::ofstream hmb("hermabin_"+interpol_str+".dat");
	for (unsigned int i=0; i< num.size(); i++) {
		hmb << int(i)-3 << '\t' << num[i] << std::endl;
	}
	std::ofstream tm("target_"+interpol_str+".dat");
	for (unsigned int x=0; x<size(0); x++) {
		for (unsigned int y=0; y<size(1); y++) {
			tm << ti[y*size(0)+x] - 1.0 << '\t';
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
	/*
		OCLH::queue.enqueueNDRangeKernel (
				rotateImg,
				cl::NullRange,
				cl::NDRange(size(0),size(1)));
				*/
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
		_data1D_rotated[i] = 0.0;
		for (hi h: _heritage_map1D[i]) {
			_data1D_rotated[i] += _data1D[h.index]*h.weight;
		}
	}
	std::copy_n(_data1D_rotated,size(0)*size(1),_data1D);
#endif // FR_USE_CL
}

void vfps::PhaseSpace::kick(const std::vector<meshdata_t>& AF)
{
	;
}

void vfps::PhaseSpace::rotateAndKick(const interpol_t deltat,
									 const std::vector<interpol_t>& AF)
{
	constexpr double o3 = 1./3.;

	const interpol_t cos_dt = std::cos(deltat);
	const interpol_t sin_dt = std::sin(deltat);

	for(unsigned int p_i=0; p_i< size(1); p_i++) {
		for (unsigned int q_i=0; q_i< size(0); q_i++) {
			//  Find cell of inverse image (qp,pp) of grid point i,j.
			interpol_t qp = cos_dt*x(0,q_i) - sin_dt*(x(1,p_i)+AF[q_i]); //q', backward mapping
			interpol_t pp = sin_dt*x(0,q_i) + cos_dt*(x(1,p_i)+AF[q_i]); //p'
			//if (willStayInMesh(qp,pp)) {
			if (insideMesh(qp,pp)) {
				unsigned int id = floor((qp-getMin(0))/getDelta(0)); //meshpoint smaller q'
				unsigned int jd = floor((pp-getMin(1))/getDelta(1)); //numper of lower mesh point from p'

				// arrays of Lagrange interpolation
				std::array<interpol_t,is> laq;
				std::array<interpol_t,is> lap;

				// gridpoint matrix used for interpolation
				std::array<std::array<meshdata_t,is>,is> ph;

				//Scaled arguments of interpolation functions:
				interpol_t xiq = (qp-x(0,id))/getDelta(0);  //distance from id
				interpol_t xip = (pp-x(1,jd))/getDelta(1); //distance of p' from lower mesh point

				for (unsigned int j1=0; j1<is; j1++) {
					unsigned int j0 = jd+j1-1;
					for (unsigned int i1=0; i1<is; i1++) {
						unsigned int i0 = id+i1-1;
						if(i0< size(0) && j0 < size(1) ){
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
						_data_rotated[q_i][p_i] += std::round(ph[i1][j1] * laq[i1] * lap[j1]);
					}
				}
				_data_rotated[q_i][p_i] /= 4;
			} else {
				_data_rotated[q_i][p_i] = 0.;
			}
		}
	}
	swapDataTmp();
}

void vfps::PhaseSpace::swapDataTmp() {
	meshdata_t** data_swap = _data;
	_data = _data_rotated;
	_data_rotated = data_swap;
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
									 sizeof(hi)*is*is*size(0)*size(1),
									 _heritage_map1D);
	applyHM2D = cl::Kernel(CLProgApplyHM::p, "applyHM2D");
	applyHM2D.setArg(0, _data_buf);
	applyHM2D.setArg(1, _heritage_map1D_buf);
	applyHM2D.setArg(2, size(1));
	applyHM2D.setArg(3, _data_rotated_buf);
	rotateImg = cl::Kernel(CLProgRotateKick::p, "rotateKick");
	rotateImg.setArg(0, _data_rotated_buf);
	rotateImg.setArg(1, _data_buf);
	rotateImg.setArg(2, img_size);
	rotateImg.setArg(3, rotation);
}
