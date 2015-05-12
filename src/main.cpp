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

	// load pattern to start with
	png::image<png::gray_pixel_16> image;
	try {
		image.read(opts.getStartDistPNG());
	} catch ( const png::std_error &e ) {
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch ( const png::error &e ) {
		std::cerr << "Problem loading " << opts.getStartDistPNG()
				  << ": " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	uint16_t ps_size = image.get_height();
	if (image.get_width() != ps_size) {
		std::cerr << "Phase space has to be quadratic. Please adjust "
				  << opts.getStartDistPNG() << std::endl;
		return EXIT_FAILURE;
	}

	PhaseSpace mesh(ps_size,-10.0,10.0,-10.0,10.0);

	for (unsigned int x=0; x<ps_size; x++) {
		for (unsigned int y=0; y<ps_size; y++) {
			mesh[x][y] = image[ps_size-y-1][x]/float(UINT16_MAX);
		}
	}

	HDF5File file(opts.getOutFile(),ps_size);

	PhaseSpace mesh_rotated(mesh);

	#ifdef INOVESA_USE_GUI
	Display display;
	#endif

	const unsigned int steps = opts.getSteps();
	const unsigned int outstep = opts.getOutSteps();
	const float rotations = opts.getNRotations();
	const double f_s = opts.getSyncFreq();
	const double t_d = opts.getDampingTime();

	/* angle of one rotation step (in rad)
	 * (angle = 2*pi corresponds to 1 synchrotron period)
	 */
	const double angle = 2*M_PI/steps;

	const double e0 = 2.0/(f_s*t_d*steps);

#define FP
#define RT
#ifdef FP
	FokkerPlanckMap fpm(&mesh_rotated,&mesh,ps_size,ps_size,
						vfps::FokkerPlanckMap::FPType::full,e0);
#else
	Identity fpm(&mesh_rotated,&mesh,ps_size,ps_size);
#endif
#ifdef RT
	RotationMap rm(&mesh,&mesh_rotated,ps_size,ps_size,angle);
#else
	Identity rm(&mesh,&mesh_rotated,ps_size,ps_size);
#endif

	#ifdef INOVESA_USE_CL
	mesh.syncCLMem(vfps::PhaseSpace::clCopyDirection::cpu2dev);
	#endif // INOVESA_USE_CL
	for (unsigned int i=0;i<steps*rotations;i++) {
		if (i%outstep == 0) {
			#ifdef INOVESA_USE_CL
			mesh.syncCLMem(vfps::PhaseSpace::clCopyDirection::dev2cpu);
			#endif // INOVESA_USE_CL
			file.append(&mesh);
			#ifdef INOVESA_USE_GUI
			display.createTexture(&mesh);
			display.draw();
			display.delTexture();
			#else
			std::cout << static_cast<float>(i)/steps << '/'
					  << rotations << std::endl;
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

