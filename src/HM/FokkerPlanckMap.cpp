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
									   const unsigned int xsize,
									   const unsigned int ysize,
									   FPType fpt, double e0) :
	  HeritageMap(in, out, xsize, ysize, 3)
{
	// the following doubles should be interpol_t
	const double e02d = e0/(2.*in->getDelta(1));
	const double e02d2 = e0/(in->getDelta(1)*in->getDelta(1));
	const double daome = 1+e0;
	const double diome = 1-2*e02d2;
	const double fptme = e0-2*e02d2+2;


	if (fpt == FPType::none) {
		for (unsigned int i=0; i< _size; i++) {
			_hinfo[i] = {0,0};
		}
	} else {
		for (unsigned int i=0; i< _xsize; i++) {
			_heritage_map[i][0][0] = {0,0};
			_heritage_map[i][0][1] = {0,0};
			_heritage_map[i][0][2] = {0,0};
			switch (fpt) {
			case FPType::damping_only:
				for (unsigned int j=1; j< _ysize-1; j++) {
					_heritage_map[i][j][0]={i*_ysize+j-1, e02d*in->x(1,j)};
					_heritage_map[i][j][1]={i*_ysize+j  , daome};
					_heritage_map[i][j][2]={i*_ysize+j+1,-e02d*in->x(1,j)};
				}
				break;
			case FPType::diffusion_only:
				for (unsigned int j=1; j< _ysize-1; j++) {
					_heritage_map[i][j][0]={i*_ysize+j-1,e02d2};
					_heritage_map[i][j][1]={i*_ysize+j  ,diome};
					_heritage_map[i][j][2]={i*_ysize+j+1,e02d2};
				}
				break;
			case FPType::full:
				for (unsigned int j=1; j< _ysize-1; j++) {
					_heritage_map[i][j][0]={i*_ysize+j-1,e02d2+e02d*in->x(1,j)};
					_heritage_map[i][j][1]={i*_ysize+j  ,fptme};
					_heritage_map[i][j][2]={i*_ysize+j+1,e02d2-e02d*in->x(1,j)};
				}
				break;
			default:
				break;
			}
			_heritage_map[i][_ysize-1][0] = {0,0};
			_heritage_map[i][_ysize-1][1] = {0,0};
			_heritage_map[i][_ysize-1][2] = {0,0};
		}
	}
}

