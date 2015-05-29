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
	HeritageMap(in,out,xsize,ysize,
				size_t(xsize)*size_t(ysize)*it*it/ROTMAP_SIZE,it*it,it),
	_sat(interpol_saturating)
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

	const meshaxis_t cos_dt = cos(angle);
	const meshaxis_t sin_dt = -sin(angle);

	#if ROTMAP_SIZE == 1
	for (meshindex_t q_i=0; q_i< _xsize; q_i++) {
	#else // ROTMAP_SIZE == 2/4
	for (meshindex_t q_i=0; q_i< _xsize/2; q_i++) {
	#endif
		for(meshindex_t p_i=0; p_i< _ysize; p_i++) {
			// Cell of inverse image (qp,pp) of grid point i,j.
			meshaxis_t qp; //q', backward mapping
			meshaxis_t pp; //p'
			// interpolation type specific q and p coordinates
			meshaxis_t pcoord=0;
			meshaxis_t qcoord=0;
			meshaxis_t qq_int;
			meshaxis_t qp_int;
			//Scaled arguments of interpolation functions:
			meshindex_t id; //meshpoint smaller q'
			meshindex_t jd; //numper of lower mesh point from p'
			interpol_t xiq; //distance from id
			interpol_t xip; //distance of p' from lower mesh point
			switch (rt) {
			case RotationCoordinates::mesh:
				qp = cos_dt*(q_i-(_xsize-1)/2.0)
						- sin_dt*(p_i-(_ysize-1)/2.0)+(_xsize-1)/2.0;
				pp = sin_dt*(q_i-(_xsize-1)/2.0)
						+ cos_dt*(p_i-(_ysize-1)/2.0)+(_ysize-1)/2.0;
				qcoord = qp;
				pcoord = pp;
				break;
			case RotationCoordinates::norm_0_1:
				qp = cos_dt*((q_i-(_xsize-1)/2.0)/(_xsize-1))
						- sin_dt*((p_i-(_ysize-1)/2.0)/(_ysize-1));
				pp = sin_dt*((q_i-(_xsize-1)/2.0)/(_xsize-1))
						+ cos_dt*((p_i-(_ysize-1)/2.0)/(_ysize-1));
				qcoord = (qp+0.5)*(_xsize-1);
				pcoord = (pp+0.5)*(_ysize-1);
				break;
			case RotationCoordinates::norm_pm1:
				qp = cos_dt*(2*static_cast<int>(q_i)-static_cast<int>(_xsize-1))
							/static_cast<meshaxis_t>(_xsize-1)
				   - sin_dt*(2*static_cast<int>(p_i)-static_cast<int>(_ysize-1))
							/static_cast<meshaxis_t>(_ysize-1);

				pp = sin_dt*(2*static_cast<int>(q_i)-static_cast<int>(_xsize-1))
							/static_cast<meshaxis_t>(_xsize-1)
				   + cos_dt*(2*static_cast<int>(p_i)-static_cast<int>(_ysize-1))
							/static_cast<meshaxis_t>(_ysize-1);
				qcoord = (qp+1)*(_xsize-1)/2;
				pcoord = (pp+1)*(_ysize-1)/2;
				break;
			}
			xiq = std::modf(qcoord, &qq_int);
			xip = std::modf(pcoord, &qp_int);
			id = qq_int;
			jd = qp_int;

			if (id <  _xsize && jd < _ysize)
			{
				// create vectors containing interpolation coefficiants
				switch(_it) {
				case InterpolationType::none:
					icq[0] = 1;

					icp[0] = 1;
					break;
				case InterpolationType::linear:
					icq[0] = interpol_t(1)-xiq;
					icq[1] = xiq;

					icp[0] = interpol_t(1)-xip;
					icp[1] = xip;
					break;
				case InterpolationType::quadratic:
					icq[0] = xiq*(xiq-interpol_t(1))/interpol_t(2);
					icq[1] = interpol_t(1)-xiq*xiq;
					icq[2] = xiq*(xiq+interpol_t(1))/interpol_t(2);

					icp[0] = xip*(xip-interpol_t(1))/interpol_t(2);
					icp[1] = interpol_t(1)-xip*xip;
					icp[2] = xip*(xip+interpol_t(1))/interpol_t(2);
					break;
				case InterpolationType::cubic:
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
				for (size_t hmq=0; hmq<_it; hmq++) {
					for (size_t hmp=0; hmp<_it; hmp++){
						hmc[hmp*_it+hmq] = icq[hmp]*icp[hmq];
					}
				}


				// renormlize to minimize rounding errors
//				renormalize(hmc.size(),hmc.data());

				// write heritage map
				for (meshindex_t j1=0; j1<_it; j1++) {
					 meshindex_t j0 = jd+j1-(_it-1)/2;
					for (meshindex_t i1=0; i1<_it; i1++) {
						 meshindex_t i0 = id+i1-(_it-1)/2;
						if(i0< _xsize && j0 < _ysize ){
							ph[i1][j1].index = i0*_ysize+j0;
							ph[i1][j1].weight = hmc[i1*_it+j1];
						} else {
							ph[i1][j1] = {0,0};
						}
						_hinfo[(q_i*ysize+p_i)*_ip+i1*_it+j1] = ph[i1][j1];
					}
				}
			}
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
				applyHM = cl::Kernel(CLProgApplyHM::p, "applyHM4sat");
				applyHM.setArg(0, _in->data_buf);
				applyHM.setArg(1, _hi_buf);
				applyHM.setArg(2, _out->data_buf);
				#elif ROTMAP_SIZE == 2
				applyHM = cl::Kernel(CLProgApplyHM::p, "applyHM4_2sat");
				applyHM.setArg(0, _in->data_buf);
				applyHM.setArg(1, _hi_buf);
				applyHM.setArg(2, _size);
				applyHM.setArg(3, _out->data_buf);
				#endif
			}
		} else {
			applyHM = cl::Kernel(CLProgApplyHM::p, "applyHM1D");
			applyHM.setArg(0, _in->data_buf);
			applyHM.setArg(1, _hi_buf);
			applyHM.setArg(2, it*it);
			applyHM.setArg(3, _out->data_buf);
		}
	}
	#endif // INOVESA_USE_CL

	delete [] icp;
	delete [] icq;
	delete [] hmc;

	delete [] ph;
	delete [] ph1D;
}

void vfps::RotationMap::apply()
{
	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		#ifdef INOVESA_SYNC_CL
		_in->syncCLMem(PhaseSpace::clCopyDirection::cpu2dev);
		#endif // INOVESA_SYNC_CL
		OCLH::queue.enqueueNDRangeKernel (
					applyHM,
					cl::NullRange,
					cl::NDRange(_size/ROTMAP_SIZE));
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
	}
}