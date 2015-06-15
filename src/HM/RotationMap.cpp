/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2014-2015: Patrik Schönfeldt                                 */
/*                                                                            */
/* This file is part of Inovesa.                                              */
/* Inovesa is free software: you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* Inovesa is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           */
/******************************************************************************/

#include "HM/RotationMap.hpp"

vfps::RotationMap::RotationMap(PhaseSpace* in, PhaseSpace* out,
							   const meshindex_t xsize,
							   const meshindex_t ysize,
							   const meshaxis_t angle,
							   const InterpolationType it,
							   const RotationCoordinates rt,
							   bool interpol_saturating) :
#if ROTMAP_SIZE == 0 // no heritage map is saved, just hi for one mesh cell
	HeritageMap(in,out,xsize,ysize,1,it*it,it),
#else
	HeritageMap(in,out,xsize,ysize,
				size_t(xsize)*size_t(ysize)*it*it/ROTMAP_SIZE,it*it,it),
#endif
	_sat(interpol_saturating),
	_rt(rt),
	_cos_dt(cos(-angle)),
	_sin_dt(sin(-angle))
{
#if ROTMAP_SIZE == 0
	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
	rot = {{_cos_dt,_sin_dt}};
	imgsize = {{cl_int(_xsize),cl_int(_ysize)}};

	genCode4Rotation();
	_cl_prog  = OCLH::prepareCLProg(_cl_code);

	applyHM = cl::Kernel(_cl_prog, "applyRotation");
	applyHM.setArg(0, _in->data_buf);
	applyHM.setArg(1, imgsize);
	applyHM.setArg(2, rot);
	applyHM.setArg(3, _out->data_buf);
	}
	#endif
#else
	#if ROTMAP_SIZE == 1
	for (meshindex_t q_i=0; q_i< _xsize; q_i++) {
	#else // ROTMAP_SIZE == 2/4
	for (meshindex_t q_i=0; q_i< _xsize/2; q_i++) {
	#endif
		for(meshindex_t p_i=0; p_i< _ysize; p_i++) {
			genHInfo(q_i,p_i,&_hinfo[(q_i*ysize+p_i)*_ip]);
		}
	}
	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		_hi_buf = cl::Buffer(OCLH::context,
							 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
							 sizeof(hi)*_ip*_size/ROTMAP_SIZE,
							 _hinfo);
		if (_sat) {
			if (it == InterpolationType::cubic) {
				#if ROTMAP_SIZE == 1
				genCode4HM4sat();
				_cl_prog  = OCLH::prepareCLProg(_cl_code);

				applyHM = cl::Kernel(_cl_prog, "applyHM4sat");
				applyHM.setArg(0, _in->data_buf);
				applyHM.setArg(1, _hi_buf);
				applyHM.setArg(2, _out->data_buf);
				#elif ROTMAP_SIZE == 2

				genCode4HM4_2sat();
				_cl_prog  = OCLH::prepareCLProg(_cl_code);

				applyHM = cl::Kernel(_cl_prog, "applyHM4_2sat");
				applyHM.setArg(0, _in->data_buf);
				applyHM.setArg(1, _hi_buf);
				applyHM.setArg(2, _size);
				applyHM.setArg(3, _out->data_buf);
				#endif
			}
		} else {
			genCode4HM1D();
			_cl_prog  = OCLH::prepareCLProg(_cl_code);

			applyHM = cl::Kernel(_cl_prog, "applyHM1D");
			applyHM.setArg(0, _in->data_buf);
			applyHM.setArg(1, _hi_buf);
			applyHM.setArg(2, _size);
			applyHM.setArg(3, _out->data_buf);
		}
	}
	#endif // INOVESA_USE_CL
#endif // ROTMAP_SIZE > 0
}

void vfps::RotationMap::apply()
{
	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		#ifdef INOVESA_SYNC_CL
		_in->syncCLMem(PhaseSpace::clCopyDirection::cpu2dev);
		#endif // INOVESA_SYNC_CL
		#if ROTMAP_SIZE == 0
		 // stay away from mesh borders
		OCLH::queue.enqueueNDRangeKernel (
					applyHM,
					cl::NDRange(1,1),
					cl::NDRange(_xsize-_it+1,_ysize-_it+1));
		#else // ROTMAP_SIZE > 0
		OCLH::queue.enqueueNDRangeKernel (
					applyHM,
					cl::NullRange,
					cl::NDRange(_size/ROTMAP_SIZE));
		#endif // ROTMAP_SIZE > 0
		#ifdef CL_VERSION_1_2
		OCLH::queue.enqueueBarrierWithWaitList();
		#else // CL_VERSION_1_2
		OCLH::queue.enqueueBarrier();
		#endif // CL_VERSION_1_2
		#ifdef INOVESA_SYNC_CL
		_out->syncCLMem(PhaseSpace::clCopyDirection::dev2cpu);
		#endif // INOVESA_SYNC_CL
	} else
	#endif // INOVESA_USE_CL
	{
		meshdata_t* data_in = _in->getData();
		meshdata_t* data_out = _out->getData();

		#if ROTMAP_SIZE == 0
			for (meshindex_t q_i=0; q_i< _xsize/2; q_i++) {
				for(meshindex_t p_i=0; p_i< _ysize; p_i++) {
					meshindex_t i = q_i*_ysize+p_i;
					data_out[i] = 0;
					data_out[_size-1-i] = 0;
					genHInfo(q_i,p_i,_hinfo);
					for (uint_fast8_t j=0; j<_ip; j++) {
						hi h = _hinfo[j];
						data_out[i] += data_in[h.index]*static_cast<meshdata_t>(h.weight);
						data_out[_size-1-i] += data_in[_size-1-h.index]*static_cast<meshdata_t>(h.weight);
					}
				}
			}
		#else // ROTMAP_SIZE > 0
		for (meshindex_t i=0; i< _size/ROTMAP_SIZE; i++) {
			data_out[i] = 0;
			#if ROTMAP_SIZE > 1
			data_out[_size-1-i] = 0;
			#endif
			for (uint_fast8_t j=0; j<_ip; j++) {
				hi h = _hinfo[i*_ip+j];
				data_out[i] += data_in[h.index]*static_cast<meshdata_t>(h.weight);
				#if ROTMAP_SIZE > 1
				data_out[_size-1-i] += data_in[_size-1-h.index]*static_cast<meshdata_t>(h.weight);
				#endif
			}
			if (_sat) {
				// handle overshooting
				meshdata_t ceil=std::numeric_limits<meshdata_t>::min();
				meshdata_t flor=std::numeric_limits<meshdata_t>::max();
				for (size_t x=1; x<=2; x++) {
					for (size_t y=1; y<=2; y++) {
						ceil = std::max(ceil,data_in[_hinfo[i*_ip+x*_it+y].index]);
						flor = std::min(flor,data_in[_hinfo[i*_ip+x*_it+y].index]);
					}
				}
				data_out[i] = std::max(std::min(ceil,data_out[i]),flor);

				#if ROTMAP_SIZE > 1
				ceil=std::numeric_limits<meshdata_t>::min();
				flor=std::numeric_limits<meshdata_t>::max();
				for (size_t x=1; x<=2; x++) {
					for (size_t y=1; y<=2; y++) {
						ceil = std::max(ceil,data_in[_size-1-_hinfo[i*_ip+x*_it+y].index]);
						flor = std::min(flor,data_in[_size-1-_hinfo[i*_ip+x*_it+y].index]);
					}
				}
				data_out[_size-1-i] = std::max(std::min(ceil,data_out[_size-1-i]),flor);
				#endif
			}
		}
		#endif
	}
}

void vfps::RotationMap::genHInfo(vfps::meshindex_t q_i,
								 vfps::meshindex_t p_i,
								 vfps::HeritageMap::hi* myhinfo)
{
	// gridpoint matrix used for interpolation
	hi* ph1D = new hi[_ip];
	hi** ph = new hi*[_it];
	for (uint_fast8_t i=0; i<_it;i++) {
		ph[i] = &ph1D[i*_it];
	}

	// arrays of interpolation coefficients
	interpol_t* icq = new interpol_t[_it];
	interpol_t* icp = new interpol_t[_it];

	interpol_t* hmc = new interpol_t[_ip];

	// Cell of inverse image (qp,pp) of grid point i,j.
	meshaxis_t qp; //q', backward mapping
	meshaxis_t pp; //p'
	// interpolation type specific q and p coordinates
	meshaxis_t pcoord;
	meshaxis_t qcoord;
	meshaxis_t qq_int;
	meshaxis_t qp_int;
	//Scaled arguments of interpolation functions:
	meshindex_t xi; //meshpoint smaller q'
	meshindex_t yi; //numper of lower mesh point from p'
	interpol_t xf; //distance from id
	interpol_t yf; //distance of p' from lower mesh point
	switch (_rt) {
	case RotationCoordinates::mesh:
		qp = _cos_dt*(q_i-(_xsize-1)/2.0)
				- _sin_dt*(p_i-(_ysize-1)/2.0)+(_xsize-1)/2.0;
		pp = _sin_dt*(q_i-(_xsize-1)/2.0)
				+ _cos_dt*(p_i-(_ysize-1)/2.0)+(_ysize-1)/2.0;
		qcoord = qp;
		pcoord = pp;
		break;
	case RotationCoordinates::norm_0_1:
		qp = _cos_dt*((q_i-(_xsize-1)/2.0)/(_xsize-1))
				- _sin_dt*((p_i-(_ysize-1)/2.0)/(_ysize-1));
		pp = _sin_dt*((q_i-(_xsize-1)/2.0)/(_xsize-1))
				+ _cos_dt*((p_i-(_ysize-1)/2.0)/(_ysize-1));
		qcoord = (qp+0.5)*(_xsize-1);
		pcoord = (pp+0.5)*(_ysize-1);
		break;
	case RotationCoordinates::norm_pm1:
		qp = _cos_dt*(2*static_cast<int>(q_i)-static_cast<int>(_xsize-1))
					/static_cast<meshaxis_t>(_xsize-1)
		   - _sin_dt*(2*static_cast<int>(p_i)-static_cast<int>(_ysize-1))
					/static_cast<meshaxis_t>(_ysize-1);

		pp = _sin_dt*(2*static_cast<int>(q_i)-static_cast<int>(_xsize-1))
					/static_cast<meshaxis_t>(_xsize-1)
		   + _cos_dt*(2*static_cast<int>(p_i)-static_cast<int>(_ysize-1))
					/static_cast<meshaxis_t>(_ysize-1);
		qcoord = (qp+1)*(_xsize-1)/2;
		pcoord = (pp+1)*(_ysize-1)/2;
		break;
	}
	xf = std::modf(qcoord, &qq_int);
	yf = std::modf(pcoord, &qp_int);
	xi = qq_int;
	yi = qp_int;

	if (xi <  _xsize && yi < _ysize) {
		// create vectors containing interpolation coefficiants
		switch(_it) {
		case InterpolationType::none:
			icq[0] = 1;

			icp[0] = 1;
			break;
		case InterpolationType::linear:
			icq[0] = interpol_t(1)-xf;
			icq[1] = xf;

			icp[0] = interpol_t(1)-yf;
			icp[1] = yf;
			break;
		case InterpolationType::quadratic:
			icq[0] = xf*(xf-interpol_t(1))/interpol_t(2);
			icq[1] = interpol_t(1)-xf*xf;
			icq[2] = xf*(xf+interpol_t(1))/interpol_t(2);

			icp[0] = yf*(yf-interpol_t(1))/interpol_t(2);
			icp[1] = interpol_t(1)-yf*yf;
			icp[2] = yf*(yf+interpol_t(1))/interpol_t(2);
			break;
		case InterpolationType::cubic:
			icq[0] = (xf-interpol_t(1))*(xf-interpol_t(2))*xf
					* interpol_t(-1./6.);
			icq[1] = (xf+interpol_t(1))*(xf-interpol_t(1))
					* (xf-interpol_t(2)) / interpol_t(2);
			icq[2] = (interpol_t(2)-xf)*xf*(xf+interpol_t(1))
					/ interpol_t(2);
			icq[3] = xf*(xf+interpol_t(1))*(xf-interpol_t(1))
					* interpol_t(1./6.);

			icp[0] = (yf-interpol_t(1))*(yf-interpol_t(2))*yf
					* interpol_t(-1./6.);
			icp[1] = (yf+interpol_t(1))*(yf-interpol_t(1))
					* (yf-interpol_t(2)) / interpol_t(2);
			icp[2] = (interpol_t(2)-yf)*yf*(yf+interpol_t(1))
					/ interpol_t(2);
			icp[3] = yf*(yf+interpol_t(1))*(yf-interpol_t(1))
					* interpol_t(1./6.);
			break;
		}

		/*  Assemble interpolation
		 * (using size_t although _it is mush smaller,
		 * so that product won't overflow)
		 */
		for (size_t hmq=0; hmq<_it; hmq++) {
			for (size_t hmp=0; hmp<_it; hmp++){
				hmc[hmp*_it+hmq] = icq[hmp]*icp[hmq];
			}
		}


		// renormlize to minimize rounding errors
		// renormalize(hmc.size(),hmc.data());

		// write heritage map
		for (meshindex_t j1=0; j1<_it; j1++) {
			 meshindex_t j0 = yi+j1-(_it-1)/2;
			for (meshindex_t i1=0; i1<_it; i1++) {
				 meshindex_t i0 = xi+i1-(_it-1)/2;
				if(i0< _xsize && j0 < _ysize ){
					ph[i1][j1].index = i0*_ysize+j0;
					ph[i1][j1].weight = hmc[i1*_it+j1];
				} else {
					ph[i1][j1] = {0,0};
				}
				myhinfo[i1*_it+j1] = ph[i1][j1];
			}
		}
	} else {
		for (uint_fast8_t i=0; i<_ip; i++) {
			myhinfo[i] = {0,0};
		}
	}

	delete [] icp;
	delete [] icq;
	delete [] hmc;

	delete [] ph;
	delete [] ph1D;
}

void vfps::RotationMap::genCode4HM4_2sat()
{
	_cl_code += R"(
	__kernel void applyHM4_2sat(const __global data_t* src,
								const __global hi* hm,
								const uint size,
								__global data_t* dst)
	{
		const uint i = get_global_id(0);
		const uint offset = i*16;
		data_t tmp;
		data_t value;
	)";
	if (std::is_same<vfps::meshdata_t,float>::value) {
		_cl_code += R"(
		data_t ceil=0.0f;
		data_t flor=10.0f;
		)";
	}
	if (std::is_same<vfps::meshdata_t,double>::value) {
		_cl_code += R"(
		data_t ceil=0.0;
		data_t flor=10.0;
		)";
	}
	#if FXP_FRACPART < 31
	if (std::is_same<vfps::meshdata_t,vfps::fixp32>::value) {
		_cl_code += R"(
		data_t ceil=INT_MIN;
		data_t flor=INT_MAX;
		)";
	}
	#endif
	if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value) {
		_cl_code += R"(
		data_t ceil=LONG_MIN;
		data_t flor=LONG_MAX;
		)";
	}
	_cl_code += R"(
		value = mult(src[hm[offset].src],hm[offset].weight);
		value += mult(src[hm[offset+1].src],hm[offset+1].weight);
		value += mult(src[hm[offset+2].src],hm[offset+2].weight);
		value += mult(src[hm[offset+3].src],hm[offset+3].weight);
		value += mult(src[hm[offset+4].src],hm[offset+4].weight);
		tmp = src[hm[offset+5].src];
		value += mult(tmp,hm[offset+5].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		tmp = src[hm[offset+6].src];
		value += mult(tmp,hm[offset+6].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		value += mult(src[hm[offset+7].src],hm[offset+7].weight);
		value += mult(src[hm[offset+8].src],hm[offset+8].weight);
		tmp = src[hm[offset+9].src];
		value += mult(tmp,hm[offset+9].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		tmp = src[hm[offset+10].src];
		value += mult(tmp,hm[offset+10].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		value += mult(src[hm[offset+11].src],hm[offset+11].weight);
		value += mult(src[hm[offset+12].src],hm[offset+12].weight);
		value += mult(src[hm[offset+13].src],hm[offset+13].weight);
		value += mult(src[hm[offset+14].src],hm[offset+14].weight);
		value += mult(src[hm[offset+15].src],hm[offset+15].weight);
		dst[i] = clamp(value,flor,ceil);
	)";

	if (std::is_same<vfps::meshdata_t,float>::value) {
		_cl_code += R"(
		ceil=0.0f;
		flor=10.0f;
		)";
	}
	if (std::is_same<vfps::meshdata_t,double>::value) {
		_cl_code += R"(
		ceil=0.0;
		flor=10.0;
		)";
	}
	#if FXP_FRACPART < 31
	if (std::is_same<vfps::meshdata_t,vfps::fixp32>::value) {
		_cl_code += R"(
		ceil=INT_MIN;
		flor=INT_MAX;
		)";
	}
	#endif
	if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value) {
		_cl_code += R"(
		ceil=LONG_MIN;
		flor=LONG_MAX;
		)";
	}

	_cl_code += R"(
		value = mult(src[size-1-hm[offset].src],hm[offset].weight);
		value += mult(src[size-1-hm[offset+1].src],hm[offset+1].weight);
		value += mult(src[size-1-hm[offset+2].src],hm[offset+2].weight);
		value += mult(src[size-1-hm[offset+3].src],hm[offset+3].weight);
		value += mult(src[size-1-hm[offset+4].src],hm[offset+4].weight);
		tmp = src[size-1-hm[offset+5].src];
		value += mult(tmp,hm[offset+5].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		tmp = src[size-1-hm[offset+6].src];
		value += mult(tmp,hm[offset+6].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		value += mult(src[size-1-hm[offset+7].src],hm[offset+7].weight);
		value += mult(src[size-1-hm[offset+8].src],hm[offset+8].weight);
		tmp = src[size-1-hm[offset+9].src];
		value += mult(tmp,hm[offset+9].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		tmp = src[size-1-hm[offset+10].src];
		value += mult(tmp,hm[offset+10].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		value += mult(src[size-1-hm[offset+11].src],hm[offset+11].weight);
		value += mult(src[size-1-hm[offset+12].src],hm[offset+12].weight);
		value += mult(src[size-1-hm[offset+13].src],hm[offset+13].weight);
		value += mult(src[size-1-hm[offset+14].src],hm[offset+14].weight);
		value += mult(src[size-1-hm[offset+15].src],hm[offset+15].weight);
		dst[size-1-i] = clamp(value,flor,ceil);
	})";

}

void vfps::RotationMap::genCode4HM4sat()
{
	_cl_code+= R"(
	__kernel void applyHM4sat(	const __global data_t* src,
								const __global hi* hm,
								__global data_t* dst)
	{
		data_t value = 0;
		const uint i = get_global_id(0);
		const uint offset = i*16;
	)";
	if (std::is_same<vfps::meshdata_t,float>::value) {
		_cl_code += R"(
			data_t ceil=0.0f;
			data_t flor=1.0f;
		)";
	}
	if (std::is_same<vfps::meshdata_t,double>::value) {
		_cl_code += R"(
			data_t ceil=0.0;
			data_t flor=1.0;
		)";
	}
	#if FXP_FRACPART < 31
	if (std::is_same<vfps::meshdata_t,vfps::fixp32>::value) {
		_cl_code += R"(
			data_t ceil=INT_MIN;
			data_t flor=INT_MAX;
		)";
	}
	#endif
	if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value) {
		_cl_code += R"(
			data_t ceil=LONG_MIN;
			data_t flor=LONG_MAX;
		)";
	}
	_cl_code += R"(
		data_t tmp;
		value += mult(src[hm[offset].src],hm[offset].weight);
		value += mult(src[hm[offset+1].src],hm[offset+1].weight);
		value += mult(src[hm[offset+2].src],hm[offset+2].weight);
		value += mult(src[hm[offset+3].src],hm[offset+3].weight);
		value += mult(src[hm[offset+4].src],hm[offset+4].weight);
		tmp = src[hm[offset+5].src];
		value += mult(tmp,hm[offset+5].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		tmp = src[hm[offset+6].src];
		value += mult(tmp,hm[offset+6].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		value += mult(src[hm[offset+7].src],hm[offset+7].weight);
		value += mult(src[hm[offset+8].src],hm[offset+8].weight);
		tmp = src[hm[offset+9].src];
		value += mult(tmp,hm[offset+9].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		tmp = src[hm[offset+10].src];
		value += mult(tmp,hm[offset+10].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		value += mult(src[hm[offset+11].src],hm[offset+11].weight);
		value += mult(src[hm[offset+12].src],hm[offset+12].weight);
		value += mult(src[hm[offset+13].src],hm[offset+13].weight);
		value += mult(src[hm[offset+14].src],hm[offset+14].weight);
		value += mult(src[hm[offset+15].src],hm[offset+15].weight);
		dst[i] = clamp(value,flor,ceil);
})";
}

void vfps::RotationMap::genCode4Rotation()
{
	if (_sat) {
		Display::printText("\tSaturation not implemented for ROTMAP_SIZE == 0");
	}
	_cl_code += R"(
	__kernel void applyRotation(	const __global data_t* src,
									const int2 imgSize,
									const float2 rot,
									__global data_t* dst)
	{
		const int x = get_global_id(0)+1;
		const int y = get_global_id(1)+1;

		const float srcx = rot.x*(x-(imgSize.x+1)/2)-rot.y*(y-(imgSize.y+1)/2)
								+(imgSize.x+1)/2;
		const float srcy = rot.y*(x-(imgSize.x+1)/2)+rot.x*(y-(imgSize.y+1)/2)
								+(imgSize.y+1)/2;
		const int xi = floor(srcx);
		const int yi = floor(srcy);
		const float xf = srcx - xi;
		const float yf = srcy - yi;
	)";

	switch (_it) {
	case InterpolationType::quadratic:
		_cl_code += R"(
		const float3 icq = (float3)(xf*(xf-1)/2,(1-xf*xf),xf*(xf+1)/2);
		const float3 icp = (float3)(yf*(yf-1)/2,(1-yf*yf),yf*(yf+1)/2);
		dst[x*imgSize.y+y] =
				icq.s0*(
					+icp.s0*src[(xi-1)*imgSize.y+yi-1]
					+icp.s1*src[(xi-1)*imgSize.y+yi]
					+icp.s2*src[(xi-1)*imgSize.y+yi+1]
				) +
				icq.s1*(
					+icp.s0*src[(xi)*imgSize.y+yi-1]
					+icp.s1*src[(xi)*imgSize.y+yi]
					+icp.s2*src[(xi)*imgSize.y+yi+1]
				) +
				icq.s2*(
					+icp.s0*src[(xi+1)*imgSize.y+yi-1]
					+icp.s1*src[(xi+1)*imgSize.y+yi]
					+icp.s2*src[(xi+1)*imgSize.y+yi+1]
				);
		})";
		break;
	case InterpolationType::cubic:
		_cl_code += R"(
		const float4 icq = (float4)(
								(xf  )*(xf-1)*(xf-2)/(-6),
								(xf+1)*(xf-1)*(xf-2)/( 2),
								(xf+1)*(xf  )*(xf-2)/(-2),
								(xf+1)*(xf  )*(xf-1)/( 6)
							);

		const float4 icp = (float4)(
								(yf  )*(yf-1)*(yf-2)/(-6),
								(yf+1)*(yf-1)*(yf-2)/( 2),
								(yf+1)*(yf  )*(yf-2)/(-2),
								(yf+1)*(yf  )*(yf-1)/( 6)
							);

		dst[x*imgSize.y+y] =
				icq.s0*(
					+icp.s0*src[(xi-1)*imgSize.y+yi-1]
					+icp.s1*src[(xi-1)*imgSize.y+yi  ]
					+icp.s2*src[(xi-1)*imgSize.y+yi+1]
					+icp.s3*src[(xi-1)*imgSize.y+yi+2]
				) +
				icq.s1*(
					+icp.s0*src[(xi  )*imgSize.y+yi-1]
					+icp.s1*src[(xi  )*imgSize.y+yi  ]
					+icp.s2*src[(xi  )*imgSize.y+yi+1]
					+icp.s3*src[(xi  )*imgSize.y+yi+2]
				) +
				icq.s2*(
					+icp.s0*src[(xi+1)*imgSize.y+yi-1]
					+icp.s1*src[(xi+1)*imgSize.y+yi  ]
					+icp.s2*src[(xi+1)*imgSize.y+yi+1]
					+icp.s3*src[(xi+1)*imgSize.y+yi+2]
				) +
				icq.s3*(
					+icp.s0*src[(xi+2)*imgSize.y+yi-1]
					+icp.s1*src[(xi+2)*imgSize.y+yi  ]
					+icp.s2*src[(xi+2)*imgSize.y+yi+1]
					+icp.s3*src[(xi+2)*imgSize.y+yi+2]
				);
		})";
		break;
	}
}
