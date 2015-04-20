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

#ifndef IDENTITY_HPP
#define IDENTITY_HPP

#include "HeritageMap.hpp"

namespace vfps
{

class Identity : public HeritageMap
{
public:
	Identity(	PhaseSpace* in, PhaseSpace* out,
				const unsigned int xsize, const unsigned int ysize) :
		HeritageMap(in, out, xsize, ysize, 0) {}

	/**
	 * @brief apply copys data from in to out
	 */
	void apply()
	{
	#ifdef INOVESA_USE_CL
	#ifdef INOVESA_SYNC_CL
	_in->syncCLMem(PhaseSpace::clCopyDirection::cpu2dev);
	#endif // INOVESA_SYNC_CL
	OCLH::queue.enqueueCopyBuffer(
				_in->data_buf, _out->data_buf,
				0,0,sizeof(meshdata_t)*ps_xsize*ps_ysize);
	#ifdef CL_VERSION_1_2
	OCLH::queue.enqueueBarrierWithWaitList();
	#else // CL_VERSION_1_2
	OCLH::queue.enqueueBarrier();
	#endif // CL_VERSION_1_2
	#ifdef INOVESA_SYNC_CL
	_out->syncCLMem(PhaseSpace::clCopyDirection::dev2cpu);
	#endif // INOVESA_SYNC_CL
	#else // INOVESA_USE_CL
	meshdata_t* data_in = _in->getData();
	meshdata_t* data_out = _out->getData();

	std::copy_n(data_in,sizeof(meshdata_t)*ps_xsize*ps_ysize,data_out);

	#endif // INOVESA_USE_CL
	}
};

}

#endif // IDENTITY_HPP
