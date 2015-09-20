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
					const unsigned int xsize, const unsigned int ysize,
					const InterpolationType it) :
	HeritageMap(in,out,xsize,ysize,it,it),
	_force(new meshaxis_t[xsize]())
{
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
		for(unsigned int p_i=0; p_i< _ysize; p_i++) {
			// interpolation type specific q and p coordinates
			meshaxis_t pcoord;
			meshaxis_t qp_int;
			//Scaled arguments of interpolation functions:
			unsigned int jd; //numper of lower mesh point from p'
			interpol_t xip; //distance of p' from lower mesh point
			pcoord = p_i+_force[q_i];
			xip = std::modf(pcoord, &qp_int);
			jd = qp_int;

			if (jd < _ysize)
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
	//				renormalize(hmc.size(),hmc.data());

				// write heritage map
				for (unsigned int j1=0; j1<_it; j1++) {
					unsigned int j0 = jd+j1-(_it-1)/2;
					if(j0 < _ysize ) {
						ph[j1].index = q_i*_ysize+j0;
						ph[j1].weight = hmc[j1];
					} else {
						ph[j1].index = 0;
						ph[j1].weight = 0;
					}
					_hinfo[(q_i*_ysize+p_i)*_ip+j1] = ph[j1];
				}
			}
		}
	}
	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		OCLH::queue.enqueueWriteBuffer
			(_hi_buf,CL_TRUE,0,
			 sizeof(meshdata_t)*_size,_hinfo);
	}
	#endif // INOVESA_USE_CL

	delete [] ph;
	delete [] hmc;

	// call original apply method
	HeritageMap::apply();
}

void vfps::KickMap::laser(meshaxis_t amplitude,
						  meshaxis_t pulselen,
						  meshaxis_t wavelen)
{
	amplitude = amplitude*meshaxis_t(_ysize/2)/_in->getMax(1);
	meshaxis_t sinarg = meshaxis_t(2*M_PI)/(wavelen*meshaxis_t(_xsize)/_in->getMax(0));
	pulselen = pulselen*meshaxis_t(_xsize)/_in->getMax(0)/meshaxis_t(2.35);
	for(meshindex_t x=0; x<_xsize; x++) {
		_force[x] +=meshaxis_t(std::exp(-std::pow(-(int(x)-int(_xsize/2)),2)))
							 /(meshaxis_t(2)*pulselen*pulselen)
				*amplitude*meshaxis_t(std::sin(double(sinarg)*x));
	}
}
