/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2007-2009: Peter Schregle (Fixed Point Math Library)         */
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

#include <climits>
#include <iostream>
#include <sstream>

#include "Display.hpp"
#include "PhaseSpace.hpp"
#include "Share.hpp"
#include "CL/CLProgs.hpp"
#include "CL/OpenCLHandler.hpp"
#include "HM/RotationMap.hpp"
#include "IO/HDF5File.hpp"

using namespace vfps;

int main(int argc, char** argv)
{
	#ifdef INOVESA_USE_CL
	// OpenCL device can be given as command line argument
	unsigned int device = 0;
	if (argc == 2 ) {
		std::stringstream dev(argv[1]);
		dev >> device;
		device--;
	}
	prepareCLEnvironment(device);
	prepareCLProgs();
	#endif // INOVESA_USE_CL

	float pattern_max;
	if (std::is_same<meshdata_t,fixp32>::value) {
		pattern_max = 1.0;
	} else {
		pattern_max = 0.5;
	}

	PhaseSpace mesh(-10.0,10.0,-10.0,10.0);

	HDF5File file("results.h5");

	// create pattern to start with
	#if TEST_PATTERN == 1
		for (unsigned int x=ps_xsize/4; x<ps_xsize*3/4; x++) {
			for (unsigned int y=ps_ysize/4; y<ps_ysize*3/4; y++) {
				mesh[x][y] = 1.0f*pattern_max;
			}
		}
		for (unsigned int x=ps_xsize/4; x<ps_xsize*3/4; x++) {
				mesh[x][ps_ysize/2] = 0.5f*pattern_max;
		}
		for (unsigned int y=ps_ysize/4; y<ps_ysize*3/4; y++) {
				mesh[ps_xsize/2][y] = 0.5f*pattern_max;
		}
		for (unsigned int x=ps_xsize/4; x<ps_xsize*3/4; x++) {
			mesh[x][x] = 0.0f;
			mesh[x][ps_ysize-x] = 0.0f;
		}
	#elif TEST_PATTERN == 2
		for (int x=0; x<int(ps_xsize); x++) {
			for (int y=0; y<int(ps_ysize); y++) {
				mesh[x][y] = std::exp(-std::pow(x-int(ps_xsize)/2,2)/patterndim_x
									  -std::pow(y-int(ps_ysize)/2,2)/patterndim_y);
			}
		}
	#elif TEST_PATTERN == 3
		for (int y=pattern_margin; y<int(ps_ysize-pattern_margin); y++) {
			for (int x=pattern_margin; x<int(ps_xsize/2); x++) {
				mesh[x][y] = 1.0f*pattern_max;
			}
		}
		for (int x=pattern_margin; x<int(ps_xsize/2); x++) {
			mesh[x][x] = 0;
			mesh[x][ps_ysize/2] = 0;
			mesh[x][ps_ysize-x] = 0;
		}
	#elif TEST_PATTERN == 4
		for (int y=pattern_margin; y<int(ps_ysize/2); y++) {
			for (int x=pattern_margin; x<int(ps_xsize/2); x++) {
				mesh[x][y] = (float(y-pattern_margin))
							/float(ps_ysize/2-pattern_margin)*pattern_max;
			}
		}
		for (int y=ps_ysize/2; y<int(ps_ysize-pattern_margin); y++) {
			for (int x=pattern_margin; x<int(ps_xsize/2); x++) {
				mesh[x][y] = 1.0f*pattern_max;
			}
		}
		for (int x=pattern_margin; x<int(ps_xsize/2); x++) {
			mesh[x][ps_ysize-x] = 0.0f*pattern_max;
		}
		for (int x=ps_xsize/2; x<int(ps_xsize-pattern_margin); x++) {
			for (int y=pattern_margin; y<int(ps_ysize/2); y++) {
				mesh[x][y] = std::exp(-std::pow(x-int(ps_xsize)/2,2)/patterndim_x
									  -std::pow(y-int(ps_ysize)/2,2)/patterndim_y)
						*pattern_max;
			}
		}
	#endif

	PhaseSpace mesh_rotated(mesh);
	file.write(&mesh_rotated);

	#ifdef INOVESA_USE_GUI
	Display display;
	display.createTexture(&mesh_rotated);
	display.draw();
	#endif

	// angle of one rotation step (in rad)
	constexpr double angle = 2*M_PI/steps;
	RotationMap rm(&mesh,&mesh_rotated,ps_xsize,ps_ysize,angle);
	unsigned int outstep = 100;
	unsigned int i;
	for (i=0;i<steps*rotations;i++) {
        if (i%outstep != 0) {
			rm.apply();
			swap(mesh,mesh_rotated);
		} else {
			#ifdef INOVESA_USE_GUI
			display.createTexture(&mesh_rotated);
			display.draw();
			#endif // INOVESA_USE_GUI

			file.append(&mesh_rotated);
			rm.apply();
			swap(mesh,mesh_rotated);
			#ifdef INOVESA_USE_GUI
			display.delTexture();
			#endif // INOVESA_USE_GUI
		}
	}
	#ifdef INOVESA_USE_CL
	OCLH::queue.flush();
	#endif // INOVESA_USE_CL

	return EXIT_SUCCESS;
}

