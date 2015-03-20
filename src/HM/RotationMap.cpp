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
							   const unsigned int xsize,
							   const unsigned int ysize,
							   const meshaxis_t angle) :
	HeritageMap(in,out,xsize,ysize)
{
	std::vector<interpol_t> ti;
	ti.resize(_xsize*_ysize);

	const meshaxis_t cos_dt = cos(angle);
	const meshaxis_t sin_dt = -sin(angle);

	for (unsigned int q_i=0; q_i< _xsize; q_i++) {
		for(unsigned int p_i=0; p_i< _ysize; p_i++) {
			// Cell of inverse image (qp,pp) of grid point i,j.
			meshaxis_t qp; //q', backward mapping
			meshaxis_t pp; //p'
			// interpolation type specific q and p coordinates
			meshaxis_t pcoord;
			meshaxis_t qcoord;
			meshaxis_t qq_int;
			meshaxis_t qp_int;
			//Scaled arguments of interpolation functions:
			unsigned int id; //meshpoint smaller q'
			unsigned int jd; //numper of lower mesh point from p'
			interpol_t xiq; //distance from id
			interpol_t xip; //distance of p' from lower mesh point
			#if ROTATION_TYPE == 1
				qp = cos_dt*(q_i-_xsize/2.0)
						- sin_dt*(p_i-_ysize/2.0)+_xsize/2.0;
				pp = sin_dt*(q_i-_xsize/2.0)
						+ cos_dt*(p_i-_ysize/2.0)+_ysize/2.0;
				qcoord = qp;
				pcoord = pp;
			#elif ROTATION_TYPE == 2
				qp = cos_dt*((q_i-_xsize/2.0)/_xsize)
						- sin_dt*((p_i-_ysize/2.0)/_ysize);
				pp = sin_dt*((q_i-_xsize/2.0)/_xsize)
						+ cos_dt*((p_i-_ysize/2.0)/_ysize);
				qcoord = (qp+0.5)*_xsize;
				pcoord = (pp+0.5)*_ysize;
			#elif ROTATION_TYPE == 3
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
			#endif
			xiq = std::modf(qcoord, &qq_int);
			xip = std::modf(pcoord, &qp_int);
			id = qq_int;
			jd = qp_int;

			if (id <  _xsize && jd < _ysize)
			{
				// gridpoint matrix used for interpolation
				std::array<std::array<hi,INTERPOL_TYPE>,INTERPOL_TYPE> ph;

				// arrays of interpolation coefficients
				std::array<interpol_t,INTERPOL_TYPE> icq;
				std::array<interpol_t,INTERPOL_TYPE> icp;

				std::array<interpol_t,INTERPOL_TYPE*INTERPOL_TYPE> hmc;

				// create vectors containing interpolation coefficiants
				#if INTERPOL_TYPE == 1
					icq[0] = 1;

					icp[0] = 1;
				#elif INTERPOL_TYPE == 2
					icq[0] = interpol_t(1)-xiq;
					icq[1] = xiq;

					icp[0] = interpol_t(1)-xip;
					icp[1] = xip;
				#elif INTERPOL_TYPE == 3
					icq[0] = xiq*(xiq-interpol_t(1))/interpol_t(2);
					icq[1] = interpol_t(1)-xiq*xiq;
					icq[2] = xiq*(xiq+interpol_t(1))/interpol_t(2);

					icp[0] = xip*(xip-interpol_t(1))/interpol_t(2);
					icp[1] = interpol_t(1)-xip*xip;
					icp[2] = xip*(xip+interpol_t(1))/interpol_t(2);
				#elif INTERPOL_TYPE == 4
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
				#endif
				//  Assemble interpolation
				for (size_t hmq=0; hmq<INTERPOL_TYPE; hmq++) {
					for (size_t hmp=0; hmp<INTERPOL_TYPE; hmp++){
						hmc[hmp*INTERPOL_TYPE+hmq] = icq[hmp]*icp[hmq];
					}
				}


				// renormlize to minimize rounding errors
//				renormalize(hmc.size(),hmc.data());

				// write heritage map
				for (unsigned int j1=0; j1<INTERPOL_TYPE; j1++) {
					unsigned int j0 = jd+j1-(INTERPOL_TYPE-1)/2;
					for (unsigned int i1=0; i1<INTERPOL_TYPE; i1++) {
						unsigned int i0 = id+i1-(INTERPOL_TYPE-1)/2;
						if(i0< _xsize && j0 < _ysize ){
							ph[i1][j1].index = i0*_ysize+j0;
							ph[i1][j1].weight = hmc[i1*INTERPOL_TYPE+j1];
						} else {
							ph[i1][j1].index = 0;
							ph[i1][j1].weight = 0;
						}
						_heritage_map[q_i][p_i][i1*INTERPOL_TYPE+j1]
								= ph[i1][j1];
					}
				}
			}
		}
	}
	#ifdef INOVESA_USE_CL
	__initOpenCL();
	#endif
}
