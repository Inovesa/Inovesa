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

#include <climits>
#include <iostream>
#include <png++/png.hpp>
#include <sstream>

#include "defines.hpp"
#include "Display.hpp"
#include "PhaseSpace.hpp"
#include "Share.hpp"
#include "CL/CLProgs.hpp"
#include "CL/OpenCLHandler.hpp"
#include "HM/FokkerPlanckMap.hpp"
#include "HM/Identity.hpp"
#include "HM/RotationMap.hpp"
#include "IO/HDF5File.hpp"
#include "IO/ProgramOptions.hpp"

using namespace vfps;

int main(int argc, char** argv)
{
	ProgramOptions opts;

	bool cont;

	try {
		cont = opts.parse(argc,argv);
	} catch(std::exception& e) {
		std::cerr << "error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	if (!cont) {
		return EXIT_SUCCESS;
	}

	#ifdef INOVESA_USE_CL
	cont = prepareCLEnvironment(opts.getCLDevice()-1);
	if (!cont) {
		return EXIT_SUCCESS;
	}
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

	// load pattern to start with
	png::image<png::gray_pixel_16> image;
	try {
		image.read(opts.getStartDistPNG());
	} catch ( const png::std_error &e ) {
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch ( const png::error &e ) {
		std::cerr << "Problem loading start.png: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	if (image.get_height() != ps_xsize ||
		image.get_width()  != ps_ysize ) {
		std::cerr	<< "Provided start.png has to have "
					<< ps_ysize << "x" << ps_xsize << " pixels."
					<< std::endl
					<< "Will now quit." << std::endl;
		return EXIT_FAILURE;
	}
	for (unsigned int x=0; x<ps_xsize; x++) {
		for (unsigned int y=0; y<ps_ysize; y++) {
			mesh[x][y] = pattern_max*(image[ps_ysize-y-1][x]
									  /float(UINT16_MAX));
		}
	}

	PhaseSpace mesh_rotated(mesh);

	#ifdef INOVESA_USE_GUI
	Display display;
	#endif


	/* angle of one rotation step (in rad)
	 * (angle = 2*pi corresponds to 1 synchrotron period)
	 */
	constexpr double angle = 2*M_PI/steps;

	constexpr double e0 = 2.0/(vfps::f_s*vfps::t_d*steps);


	FokkerPlanckMap fpm(&mesh_rotated,&mesh,ps_xsize,ps_ysize,
						vfps::FokkerPlanckMap::FPType::full,e0);
	RotationMap rm(&mesh,&mesh_rotated,ps_xsize,ps_ysize,angle);

	#ifdef INOVESA_USE_CL
	mesh.syncCLMem(vfps::PhaseSpace::clCopyDirection::cpu2dev);
	#endif // INOVESA_USE_CL
	unsigned int outstep = 100;
	for (unsigned int i=0;i<steps*vfps::rotations;i++) {
		if (i%outstep == 0) {
			#ifdef INOVESA_USE_CL
			mesh.syncCLMem(vfps::PhaseSpace::clCopyDirection::dev2cpu);
			#endif // INOVESA_USE_CL
			file.append(&mesh);
			#ifdef INOVESA_USE_GUI
			display.createTexture(&mesh);
			display.draw();
			display.delTexture();
			#endif // INOVESA_USE_GUI
		}
		rm.apply();
		fpm.apply();
	}
	#ifdef INOVESA_USE_CL
	OCLH::queue.flush();
	#endif // INOVESA_USE_CL

	return EXIT_SUCCESS;
}

