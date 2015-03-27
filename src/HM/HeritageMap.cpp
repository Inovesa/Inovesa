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

#include "HM/HeritageMap.hpp"

vfps::HeritageMap::HeritageMap(PhaseSpace* in, PhaseSpace* out,
							   unsigned int xsize, unsigned int ysize,
							   unsigned int interpoints) :
	_ip(interpoints),
	_heritage_map(new hi**[xsize]),
	_heritage_map1D(new hi*[xsize*ysize]()),
	_hinfo(new hi[xsize*ysize*interpoints]),
	_size(xsize*ysize),
	_xsize(xsize),
	_ysize(ysize),
	_in(in),
	_out(out)
{
	for (unsigned int i=0; i<xsize; i++) {
		for (unsigned int j=0; j<ysize; j++) {
			_heritage_map1D[i*ysize+j]=&(_hinfo[(i*ysize+j)*interpoints]);
		}
		_heritage_map[i] = &(_heritage_map1D[i*ysize]);
	}
}

vfps::HeritageMap::~HeritageMap()
{
	delete [] _heritage_map1D;
	delete [] _heritage_map;
	delete [] _hinfo;
}

void vfps::HeritageMap::apply()
{
	#ifdef INOVESA_USE_CL
	OCLH::queue.enqueueWriteBuffer
				(_in->data_buf, CL_TRUE,
				 0,sizeof(meshdata_t)*ps_xsize*ps_ysize,
				_in->getData());
	OCLH::queue.enqueueNDRangeKernel (
				applyHM,
				cl::NullRange,
				cl::NDRange(_size));
	#ifdef CL_VERSION_1_2
	OCLH::queue.enqueueBarrierWithWaitList();
	#else // CL_VERSION_1_2
	OCLH::queue.enqueueBarrier();
	#endif // CL_VERSION_1_2
	OCLH::queue.enqueueReadBuffer
				(_out->data_buf, CL_TRUE,
				 0,sizeof(meshdata_t)*_size,
				_out->getData());
	#else // INOVESA_USE_CL
	meshdata_t* data_in = _in->getData();
	meshdata_t* data_out = _out->getData();

	for (unsigned int i=0; i< _size; i++) {
		data_out[i] = 0;
		for (unsigned int j=0; j<_ip; j++) {
			hi h = _heritage_map1D[i][j];
			data_out[i] += data_in[h.index]*static_cast<meshdata_t>(h.weight);
		}
	}
	#endif // INOVESA_USE_CL
}

#ifdef INOVESA_USE_CL
void vfps::HeritageMap::__initOpenCL()
{
	_hi_buf = cl::Buffer(OCLH::context,
									 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
									 sizeof(hi)*_ip*_size,
									 _hinfo);
	#if INTERPOL_TYPE == 4
		applyHM = cl::Kernel(CLProgApplyHM::p, "applyHM4sat");
		applyHM.setArg(0, _in->data_buf);
		applyHM.setArg(1, _hi_buf);
		applyHM.setArg(2, _out->data_buf);
	#else
		applyHM = cl::Kernel(CLProgApplyHM::p, "applyHM1D");
		applyHM.setArg(0, _in->data_buf);
		applyHM.setArg(1, _heritage_map1D_buf);
		applyHM.setArg(2, INTERPOL_TYPE*INTERPOL_TYPE);
		applyHM.setArg(3, _out->data_buf);
	#endif
}
#endif
