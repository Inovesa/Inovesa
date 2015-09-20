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

#include "HM/WakeKickMap.hpp"

vfps::WakeKickMap::WakeKickMap(vfps::PhaseSpace* in, vfps::PhaseSpace* out,
				const vfps::meshindex_t xsize, const vfps::meshindex_t ysize,
				const std::pair<vfps::meshindex_t,double>* wakefunction,
				const size_t wakesize, const InterpolationType it) :
	KickMap(in,out,xsize,ysize,it),
	_wakefunction(new meshaxis_t[2*xsize])
{
	for (size_t i=0; i<wakesize && i < 2*xsize; i++) {
		_wakefunction[i] = wakefunction[i].second;
	}
	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		genCode4HM1D();
		_cl_prog  = OCLH::prepareCLProg(_cl_code);

		_hi_buf = cl::Buffer(OCLH::context,
							 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
							 sizeof(hi)*_ip*_size,
							 _hinfo);

		applyHM = cl::Kernel(_cl_prog, "applyHM1D");
		applyHM.setArg(0, _in->data_buf);
		applyHM.setArg(1, _hi_buf);
		applyHM.setArg(2, _size);
		applyHM.setArg(3, _out->data_buf);
	}
	#endif // INOVESA_USE_CL
}

vfps::WakeKickMap::~WakeKickMap()
{
	delete [] _wakefunction;
}

void vfps::WakeKickMap::apply()
{
	integral_t charge = _in->integral();
	const integral_t* density = _in->getProjection(0);
	for (unsigned int i=0;i<_xsize;i++) {
		_force[i] = 0;
		for (unsigned int j=0;j<_xsize;j++) {
			_force[i] += meshaxis_t(density[j]/charge*_wakefunction[_xsize+i-j]);
		}
	}

	KickMap::apply();
}

