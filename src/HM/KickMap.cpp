/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2015: Patrik Schönfeldt                                 *
 *                                                                            *
 * This file is part of Inovesa.                                              *
 * Inovesa is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Inovesa is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           *
 ******************************************************************************/

#include "HM/KickMap.hpp"

vfps::KickMap::KickMap(vfps::PhaseSpace* in, vfps::PhaseSpace* out,
					const meshindex_t xsize, const meshindex_t ysize,
					const InterpolationType it) :
	HeritageMap(in,out,xsize,1,it,it),
	_force(new meshaxis_t[xsize]()),
	_meshysize(ysize)
{
	#ifdef INOVESA_INIT_KICKMAP
	for (meshindex_t q_i=0; q_i<xsize; q_i++) {
		_hinfo[q_i*_ip].index = _meshysize/2;
		_hinfo[q_i*_ip].weight = 1;
		for (unsigned int j1=1; j1<_it; j1++) {
			_hinfo[q_i*_ip+j1].index = 0;
			_hinfo[q_i*_ip+j1].weight = 0;
		}
	}
	#endif
	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		_cl_code += R"(
		__kernel void applyHM_Kick(	const __global data_t* src,
									const __global hi* hm,
									const uint hm_len,
									const uint ysize,
									__global data_t* dst)
		{
			data_t value = 0;
			const uint x = get_global_id(0);
			const uint y = get_global_id(1);
			const uint hmoffset = x*hm_len;
			const uint meshoffs = x*ysize;
			for (uint j=0; j<hm_len; j++)
			{
				value += mult(	src[meshoffs+y+hm[hmoffset+j].src-ysize/2],
								hm[hmoffset+j].weight);
			}
			dst[meshoffs+y] = value;
		}
		)";
		_cl_prog  = OCLH::prepareCLProg(_cl_code);

		_hi_buf = cl::Buffer(OCLH::context,
							 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
							 sizeof(hi)*_ip*_size,
							 _hinfo);

		applyHM = cl::Kernel(_cl_prog, "applyHM_Kick");
		applyHM.setArg(0, _in->data_buf);
		applyHM.setArg(1, _hi_buf);
		applyHM.setArg(2, _ip);
		applyHM.setArg(3, _meshysize);
		applyHM.setArg(4, _out->data_buf);
	}
	#endif // INOVESA_USE_CL
}

vfps::KickMap::~KickMap()
{
	delete [] _force;
}

void vfps::KickMap::apply()
{
	// gridpoint matrix used for interpolation
	hi* ph = new hi[_it];

	// arrays of interpolation coefficients
	interpol_t* hmc = new interpol_t[_it];

	// translate force into HM
	for (unsigned int q_i=0; q_i< _xsize; q_i++) {
		meshaxis_t poffs;
		meshaxis_t qp_int;
		//Scaled arguments of interpolation functions:
		meshindex_t jd; //numper of lower mesh point from p'
		interpol_t xip; //distance of p' from lower mesh point
		poffs = _meshysize/2+_force[q_i];
		xip = std::modf(poffs, &qp_int);
		jd = qp_int;

		if (jd < _meshysize)
		{
			// create vectors containing interpolation coefficiants
			switch (_it) {
			case InterpolationType::none:
				hmc[0] = 1;
				break;
			case InterpolationType::linear:
				hmc[0] = interpol_t(1)-xip;
				hmc[1] = xip;
				break;
			case InterpolationType::quadratic:
				hmc[0] = xip*(xip-interpol_t(1))/interpol_t(2);
				hmc[1] = interpol_t(1)-xip*xip;
				hmc[2] = xip*(xip+interpol_t(1))/interpol_t(2);
				break;
			case InterpolationType::cubic:
				hmc[0] = (xip-interpol_t(1))*(xip-interpol_t(2))*xip
						* interpol_t(-1./6.);
				hmc[1] = (xip+interpol_t(1))*(xip-interpol_t(1))
						* (xip-interpol_t(2)) / interpol_t(2);
				hmc[2] = (interpol_t(2)-xip)*xip*(xip+interpol_t(1))
						/ interpol_t(2);
				hmc[3] = xip*(xip+interpol_t(1))*(xip-interpol_t(1))
						* interpol_t(1./6.);
				break;
			}

			// renormlize to minimize rounding errors
			// renormalize(hmc.size(),hmc.data());

			// write heritage map
			for (unsigned int j1=0; j1<_it; j1++) {
				unsigned int j0 = jd+j1-(_it-1)/2;
				if(j0 < _meshysize ) {
					ph[j1].index = j0;
					ph[j1].weight = hmc[j1];
				} else {
					ph[j1].index = _meshysize/2;
					ph[j1].weight = 0;
				}
				_hinfo[q_i*_ip+j1] = ph[j1];
			}
		} else {
			for (unsigned int j1=0; j1<_it; j1++) {
				ph[j1].index = _meshysize/2;
				ph[j1].weight = 0;
				_hinfo[q_i*_ip+j1] = ph[j1];
			}
		}
	}

	delete [] ph;
	delete [] hmc;

	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		OCLH::queue.enqueueWriteBuffer
			(_hi_buf,CL_TRUE,0,
			 sizeof(hi)*_ip*_size,_hinfo);
		#ifdef INOVESA_SYNC_CL
		_in->syncCLMem(PhaseSpace::clCopyDirection::cpu2dev);
		#endif // INOVESA_SYNC_CL
		OCLH::queue.enqueueNDRangeKernel (
					applyHM,
					cl::NullRange,
					cl::NDRange(_xsize,_meshysize));
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


	for (meshindex_t x=0; x< _xsize; x++) {
		const meshindex_t offs = x*_meshysize;
		for (meshindex_t y=0; y< _meshysize; y++) {
			data_out[offs+y] = 0;
			for (uint_fast8_t j=0; j<_ip; j++) {
				hi h = _hinfo[x*_ip+j];
				meshindex_t ysrc = static_cast<int32_t>(y+h.index)
								 - static_cast<int32_t>(_meshysize/2);
				if (ysrc < _meshysize) {
					data_out[offs+y] += data_in[offs+ysrc]
								*static_cast<meshdata_t>(h.weight);
				}
			}
		}
	}
	}
}

void vfps::KickMap::laser(meshaxis_t amplitude,
						  meshaxis_t pulselen,
						  meshaxis_t wavelen)
{
	amplitude = amplitude*meshaxis_t(_meshysize/2)/_in->getMax(1);
	meshaxis_t sinarg = meshaxis_t(2*M_PI)/(wavelen*meshaxis_t(_xsize)/_in->getMax(0));
	pulselen = pulselen*meshaxis_t(_xsize)/_in->getMax(0)/meshaxis_t(2.35);
	for(meshindex_t x=0; x<_xsize; x++) {
		_force[x] +=meshaxis_t(std::exp(-std::pow(-(int(x)-int(_xsize/2)),2)))
							 /(meshaxis_t(2)*pulselen*pulselen)
				*amplitude*meshaxis_t(std::sin(double(sinarg)*x));
	}
}
