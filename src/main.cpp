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

constexpr meshaxis_t z0=0;
constexpr meshaxis_t sz=1;
constexpr meshaxis_t sd=1;

meshdata_t gaus(meshdata_t x, meshdata_t x0, meshdata_t sx)
{
	return std::exp(-(((x-x0)*(x-x0))/(2*sx*sx)));
}
meshdata_t gaus2d(meshdata_t x, meshdata_t x0, meshdata_t sx,
				  meshdata_t y, meshdata_t y0, meshdata_t sy)
{
	return std::exp(-(((x-x0)*(x-x0))/(2*sx*sx)+((y-y0)*(y-y0)))/(2*sy*sy));
}

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

	uint32_t ps_size = 0x19000;
//	uint32_t ps_size = 0x100;

	PhaseSpace mesh(ps_size,-10.0,10.0,-10.0,10.0);

	for (unsigned int x=0; x<ps_size; x++) {
		for (unsigned int y=0; y<ps_size; y++) {
			meshaxis_t z = mesh.x(0,x);
			meshaxis_t d = mesh.x(1,y);
			meshaxis_t dd = 7*gaus(z,0,1e-3);
			mesh[x][y] = 0.5*gaus2d(z,z0,sz,d, dd,sd)
					   + 0.5*gaus2d(z,z0,sz,d,-dd,sd);
		}
	}

	HDF5File file(opts.getOutFile(),ps_size);

	PhaseSpace mesh_rotated(mesh);

	#ifdef INOVESA_USE_GUI
	Display* display = nullptr;
	if (opts.showPhaseSpace()) {
		display = new Display();
	}
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
			if (opts.showPhaseSpace()) {
				display->createTexture(&mesh);
				display->draw();
				display->delTexture();
			}
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

