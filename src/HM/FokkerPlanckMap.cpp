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

#include "HM/FokkerPlanckMap.hpp"

vfps::FokkerPlanckMap::FokkerPlanckMap(PhaseSpace* in, PhaseSpace* out,
									   const meshindex_t xsize,
									   const meshindex_t ysize,
									   FPType fpt, double e1)
	:
	#if DERIVATION_TYPE == 1 // have to use type 2 for second derivative
	  HeritageMap(in, out, xsize, ysize, 3)
	#else
	  HeritageMap(in, out, xsize, ysize, DERIVATION_TYPE+1)
	#endif
{
	#if DERIVATION_TYPE == 1
	// the following doubles should be interpol_t
	const double e1_1d = e1/(in->getDelta(1));
	const double e1_d2 = e1/(in->getDelta(1)*in->getDelta(1));
	for (meshindex_t i=0; i< _xsize; i++) {
		_heritage_map[i][0][0] = {0,0};
		_heritage_map[i][0][1] = {0,0};
		_heritage_map[i][0][2] = {0,0};
		for (meshindex_t j=1; j< _ysize-2; j++) {
			_heritage_map[i][j][0]={i*_ysize+j-1,0};
			_heritage_map[i][j][1]={i*_ysize+j  ,1};
			_heritage_map[i][j][2]={i*_ysize+j+1,0};
		}
		if (fpt == FPType::full || fpt == FPType::damping_only) {
			for (uint16_t j=1; j<_ysize/2; j++) {
				const double pos = in->x(1,j);
				_heritage_map[i][j][0].weight +=   -e1_1d*pos;
				_heritage_map[i][j][1].weight += e1+e1_1d*pos;
			}
			for (uint16_t j=_ysize/2; j< _ysize-1; j++) {
				const double pos = in->x(1,j);
				_heritage_map[i][j][1].weight += e1-e1_1d*pos;
				_heritage_map[i][j][2].weight +=    e1_1d*pos;
			}
		}
		for (meshindex_t j=1; j< _ysize-2; j++) {
			if (fpt == FPType::full || fpt == FPType::diffusion_only) {
				_heritage_map[i][j][0].weight +=    e1_d2;
				_heritage_map[i][j][1].weight += -2*e1_d2;
				_heritage_map[i][j][2].weight +=    e1_d2;
			}
		}
		_heritage_map[i][_ysize-1][0] = {0,0};
		_heritage_map[i][_ysize-1][1] = {0,0};
		_heritage_map[i][_ysize-1][2] = {0,0};
	}
	#elif DERIVATION_TYPE == 2
	// the following doubles should be interpol_t
	const double e1_2d = e1/(2.*in->getDelta(1));
	const double e1_d2 = e1/(in->getDelta(1)*in->getDelta(1));
	for (meshindex_t i=0; i< _xsize; i++) {
		_heritage_map[i][0][0] = {0,0};
		_heritage_map[i][0][1] = {0,0};
		_heritage_map[i][0][2] = {0,0};
		for (meshindex_t j=1; j< _ysize-1; j++) {
			_heritage_map[i][j][0]={i*_ysize+j-1,0};
			_heritage_map[i][j][1]={i*_ysize+j  ,1};
			_heritage_map[i][j][2]={i*_ysize+j+1,0};

			if (fpt == FPType::full || fpt == FPType::damping_only) {
				const double pos = in->x(1,j);
				_heritage_map[i][j][0].weight += -e1_2d*pos;
				_heritage_map[i][j][1].weight +=  e1;
				_heritage_map[i][j][2].weight += +e1_2d*pos;
			}
			if (fpt == FPType::full || fpt == FPType::diffusion_only) {
				_heritage_map[i][j][0].weight +=    e1_d2;
				_heritage_map[i][j][1].weight += -2*e1_d2;
				_heritage_map[i][j][2].weight +=    e1_d2;
			}
		}
		_heritage_map[i][_ysize-1][0] = {0,0};
		_heritage_map[i][_ysize-1][1] = {0,0};
		_heritage_map[i][_ysize-1][2] = {0,0};
	}
	#elif DERIVATION_TYPE == 3
	// the following doubles should be interpol_t
	const double e1_6d = e1/(6.*static_cast<double>(in->getDelta(1)));
	const double e1_d2 = e1/(in->getDelta(1)*in->getDelta(1));
	for (meshindex_t i=0; i< _xsize; i++) {
		_heritage_map[i][0][0] = {0,0};
		_heritage_map[i][0][1] = {0,0};
		_heritage_map[i][0][2] = {0,0};
		_heritage_map[i][0][3] = {0,0};
		_heritage_map[i][1][0] = {0,0};
		_heritage_map[i][1][1] = {0,0};
		_heritage_map[i][1][2] = {0,0};
		_heritage_map[i][1][3] = {0,0};
		for (meshindex_t j=2; j< _ysize/2; j++) {
			const double pos = in->x(1,j);
			_heritage_map[i][j][0]={i*_ysize+j-2,0};
			_heritage_map[i][j][1]={i*_ysize+j-1,0};
			_heritage_map[i][j][2]={i*_ysize+j  ,1};
			_heritage_map[i][j][3]={i*_ysize+j+1,0};
			if (fpt == FPType::full || fpt == FPType::damping_only) {
				_heritage_map[i][j][0].weight +=    e1_6d*( 1.)*pos;
				_heritage_map[i][j][1].weight +=    e1_6d*(-6.)*pos;
				_heritage_map[i][j][2].weight += e1+e1_6d*( 3.)*pos;
				_heritage_map[i][j][3].weight +=    e1_6d*( 2.)*pos;
			}
			if (fpt == FPType::full || fpt == FPType::diffusion_only) {
				_heritage_map[i][j][1].weight +=    e1_d2;
				_heritage_map[i][j][2].weight += -2*e1_d2;
				_heritage_map[i][j][3].weight +=    e1_d2;
			}
		}
		for (meshindex_t j=_ysize/2; j<static_cast<meshindex_t>(_ysize-2);j++) {
			const double pos = in->x(1,j);
			_heritage_map[i][j][0]={i*_ysize+j-1,0};
			_heritage_map[i][j][1]={i*_ysize+j  ,1};
			_heritage_map[i][j][2]={i*_ysize+j+1,0};
			_heritage_map[i][j][3]={i*_ysize+j+2,0};

			if (fpt == FPType::full || fpt == FPType::damping_only) {
				_heritage_map[i][j][0].weight +=    e1_6d*(-2.)*pos;
				_heritage_map[i][j][1].weight += e1+e1_6d*(-3.)*pos;
				_heritage_map[i][j][2].weight +=    e1_6d*( 6.)*pos;
				_heritage_map[i][j][3].weight +=    e1_6d*(-1.)*pos;
			}
			if (fpt == FPType::full || fpt == FPType::diffusion_only) {
				_heritage_map[i][j][0].weight +=    e1_d2;
				_heritage_map[i][j][1].weight += -2*e1_d2;
				_heritage_map[i][j][2].weight +=    e1_d2;
			}
		}
		_heritage_map[i][_ysize-2][0] = {0,0};
		_heritage_map[i][_ysize-2][1] = {0,0};
		_heritage_map[i][_ysize-2][2] = {0,0};
		_heritage_map[i][_ysize-2][3] = {0,0};
		_heritage_map[i][_ysize-1][0] = {0,0};
		_heritage_map[i][_ysize-1][1] = {0,0};
		_heritage_map[i][_ysize-1][2] = {0,0};
		_heritage_map[i][_ysize-1][3] = {0,0};
	}
	#endif

	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		_hi_buf = cl::Buffer(OCLH::context,
							 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
							 sizeof(hi)*_ip*_size,
							 _hinfo);
		applyHM = cl::Kernel(CLProgApplyHM::p, "applyHM1D");
		applyHM.setArg(0, _in->data_buf);
		applyHM.setArg(1, _hi_buf);
		applyHM.setArg(2, _ip);
		applyHM.setArg(3, _out->data_buf);
	}
	#endif
}

