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
							   size_t xsize, size_t ysize,
							   uint_fast8_t interpoints,
							   uint_fast8_t intertype) :
	_ip(interpoints),
	_it(intertype),
	_hinfo(new hi[xsize*ysize*interpoints]),
	_size(xsize*ysize),
	_xsize(xsize),
	_ysize(ysize),
	_in(in),
	_out(out)
{
}

vfps::HeritageMap::~HeritageMap()
{
	delete [] _hinfo;
}

void vfps::HeritageMap::apply()
{
	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		#ifdef INOVESA_SYNC_CL
		_in->syncCLMem(PhaseSpace::clCopyDirection::cpu2dev);
		#endif // INOVESA_SYNC_CL
		OCLH::queue.enqueueNDRangeKernel (
					applyHM,
					cl::NullRange,
					cl::NDRange(_size));
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

		for (meshindex_t i=0; i< _size; i++) {
			data_out[i] = 0;
			for (meshindex_t j=0; j<_ip; j++) {
				hi h = _hinfo[i*_ip+j];
				data_out[i] += data_in[h.index]*static_cast<meshdata_t>(h.weight);
			}
		}
	}
}
