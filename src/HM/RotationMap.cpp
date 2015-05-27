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
	HeritageMap(in,out,xsize,ysize,it*it,it),
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

	for (meshindex_t q_i=0; q_i< _xsize/2; q_i++) {
		for(meshindex_t p_i=0; p_i< _ysize; p_i++) {
			// Cell of inverse image (qp,pp) of grid point i,j.
			meshaxis_t qp; //q', backward mapping
			meshaxis_t pp; //p'
			// interpolation type specific q and p coordinates
			meshaxis_t pcoord;
			meshaxis_t qcoord;
			meshaxis_t qq_int;
			meshaxis_t qp_int;
			//Scaled arguments of interpolation functions:
			meshindex_t id; //meshpoint smaller q'
			meshindex_t jd; //numper of lower mesh point from p'
			interpol_t xiq; //distance from id
			interpol_t xip; //distance of p' from lower mesh point
			switch (rt) {
			case RotationCoordinates::mesh:
				qp = cos_dt*(q_i-_xsize/2.0)
						- sin_dt*(p_i-_ysize/2.0)+_xsize/2.0;
				pp = sin_dt*(q_i-_xsize/2.0)
						+ cos_dt*(p_i-_ysize/2.0)+_ysize/2.0;
				qcoord = qp;
				pcoord = pp;
				break;
			case RotationCoordinates::norm_0_1:
				qp = cos_dt*((q_i-_xsize/2.0)/_xsize)
						- sin_dt*((p_i-_ysize/2.0)/_ysize);
				pp = sin_dt*((q_i-_xsize/2.0)/_xsize)
						+ cos_dt*((p_i-_ysize/2.0)/_ysize);
				qcoord = (qp+0.5)*_xsize;
				pcoord = (pp+0.5)*_ysize;
				break;
			case RotationCoordinates::norm_pm1:
				qp = cos_dt*meshaxis_t(
							(2*static_cast<int>(q_i)-static_cast<int>(_xsize))
							/static_cast<int>(_xsize))
				   - sin_dt*meshaxis_t(
							(2*static_cast<int>(p_i)-static_cast<int>(_ysize))
							/static_cast<int>(_ysize));

				pp = sin_dt*meshaxis_t(
							(2*static_cast<int>(q_i)-static_cast<int>(_xsize))
							/static_cast<int>(_xsize))
				   + cos_dt*meshaxis_t(
							(2*static_cast<int>(p_i)-static_cast<int>(_ysize))
							/static_cast<int>(_ysize));
				qcoord = (qp+1)*_xsize/2;
				pcoord = (pp+1)*_ysize/2;
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
							 sizeof(hi)*_ip*_size,
							 _hinfo);
		if (_sat) {
			if (it == InterpolationType::cubic) {
				applyHM = cl::Kernel(CLProgApplyHM::p, "applyHM4sat");
				applyHM.setArg(0, _in->data_buf);
				applyHM.setArg(1, _hi_buf);
				applyHM.setArg(2, _out->data_buf);
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
		HeritageMap::apply();
	} else
	#endif // INOVESA_USE_CL
	{
		meshdata_t* data_in = _in->getData();
		meshdata_t* data_out = _out->getData();

		for (meshindex_t i=0; i< _size/2; i++) {
			data_out[i] = 0;
			data_out[_size-1-i] = 0;
			for (uint_fast8_t j=0; j<_ip; j++) {
				hi h = _hinfo[i*_ip+j];
				data_out[i] += data_in[h.index]*static_cast<meshdata_t>(h.weight);
				data_out[_size-1-i] += data_in[_size-1-h.index]*static_cast<meshdata_t>(h.weight);
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


				ceil=std::numeric_limits<meshdata_t>::min();
				flor=std::numeric_limits<meshdata_t>::max();
				for (size_t x=1; x<=2; x++) {
					for (size_t y=1; y<=2; y++) {
						ceil = std::max(ceil,data_in[_size-1-_hinfo[i*_ip+x*_it+y].index]);
						flor = std::min(flor,data_in[_size-1-_hinfo[i*_ip+x*_it+y].index]);
					}
				}
				data_out[_size-1-i] = std::max(std::min(ceil,data_out[_size-1-i]),flor);
			}
		}
	}
}
