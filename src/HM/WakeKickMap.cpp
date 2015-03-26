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
					const unsigned int xsize, const unsigned int ysize,
					const std::array<vfps::integral_t,2*ps_xsize> wake) :
	HeritageMap(in,out,xsize,ysize,INTERPOL_TYPE),
	_wake(wake)
{
}

void vfps::WakeKickMap::apply()
{
	integral_t* density = _in->projectionToX();
	for (unsigned int i=0;i<_xsize;i++) {
		_wakeforce[i] = 0;
		for (unsigned int j=0;j<_xsize;j++) {
			_wakeforce[i] += meshaxis_t(density[j]*_wake[_xsize+i-j]);
		}
	}
	// translate force into HM
	for (unsigned int q_i=0; q_i< _xsize; q_i++) {
		for(unsigned int p_i=0; p_i< _ysize; p_i++) {
			// interpolation type specific q and p coordinates
			meshaxis_t pcoord;
			meshaxis_t qp_int;
			//Scaled arguments of interpolation functions:
			unsigned int jd; //numper of lower mesh point from p'
			interpol_t xip; //distance of p' from lower mesh point
			pcoord = _wakeforce[q_i];
			xip = std::modf(pcoord, &qp_int);
			jd = qp_int;

			if (jd < _ysize)
			{
				// gridpoint matrix used for interpolation
				std::array<hi,INTERPOL_TYPE> ph;

				// arrays of interpolation coefficients
				std::array<interpol_t,INTERPOL_TYPE> hmc;

				// create vectors containing interpolation coefficiants
				#if INTERPOL_TYPE == 1
					icp[0] = 1;
				#elif INTERPOL_TYPE == 2
					icp[0] = interpol_t(1)-xip;
					icp[1] = xip;
				#elif INTERPOL_TYPE == 3
					icp[0] = xip*(xip-interpol_t(1))/interpol_t(2);
					icp[1] = interpol_t(1)-xip*xip;
					icp[2] = xip*(xip+interpol_t(1))/interpol_t(2);
				#elif INTERPOL_TYPE == 4
					hmc[0] = (xip-interpol_t(1))*(xip-interpol_t(2))*xip
							* interpol_t(-1./6.);
					hmc[1] = (xip+interpol_t(1))*(xip-interpol_t(1))
							* (xip-interpol_t(2)) / interpol_t(2);
					hmc[2] = (interpol_t(2)-xip)*xip*(xip+interpol_t(1))
							/ interpol_t(2);
					hmc[3] = xip*(xip+interpol_t(1))*(xip-interpol_t(1))
							* interpol_t(1./6.);
				#endif

				// renormlize to minimize rounding errors
//				renormalize(hmc.size(),hmc.data());

				// write heritage map
				for (unsigned int j1=0; j1<INTERPOL_TYPE; j1++) {
					unsigned int j0 = jd+j1-(INTERPOL_TYPE-1)/2;
					if(j0 < _ysize ) {
						ph[j0].index = q_i*_ysize+j0;
						ph[j0].weight = hmc[j1];
					} else {
						ph[j0].index = 0;
						ph[j0].weight = 0;
					}
					_heritage_map[q_i][p_i][j1]
							= ph[j1];
				}
			}
		}
	}

	// call original apply method
	HeritageMap::apply();
}

