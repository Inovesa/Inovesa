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

#include <chrono>
#include <climits>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#ifdef INOVESA_USE_PNG
#include <png++/png.hpp>
#endif
#include <random>
#include <sstream>

#include "defines.hpp"
#include "IO/Display.hpp"
#include "PhaseSpace.hpp"
#include "CL/OpenCLHandler.hpp"
#include "HM/FokkerPlanckMap.hpp"
#include "HM/Identity.hpp"
#include "HM/KickMap.hpp"
#include "HM/RotationMap.hpp"
#include "IO/HDF5File.hpp"
#include "IO/ProgramOptions.hpp"

using namespace vfps;

int main(int argc, char** argv)
{
	Display::start_time = std::chrono::system_clock::now();

	std::time_t start_ctime
			= std::chrono::system_clock::to_time_t(Display::start_time);
	std::stringstream sstream;
	sstream << std::ctime(&start_ctime);
	std::string timestring = sstream.str();
	timestring.resize(timestring.size()-1);

	sstream.str("");
	sstream << INOVESA_VERSION_RELEASE << '.'
			<< INOVESA_VERSION_MINOR << '.'
			<< INOVESA_VERSION_FIX;

	Display::printText("Started Inovesa (v"
					   +sstream.str()+") at "+timestring);

	ProgramOptions opts;

	try {
		if (!opts.parse(argc,argv)) {
			return EXIT_SUCCESS;
		}
	} catch(std::exception& e) {
		std::cerr << "error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	#ifdef INOVESA_USE_CL
	OCLH::active = (opts.getCLDevice() >= 0);
	if (OCLH::active) {
		try {
			OCLH::prepareCLEnvironment();
		} catch (cl::Error& e) {
			Display::printText(e.what());
			Display::printText("Will fall back to sequential version.");
			OCLH::active = false;
		}
	}
	if (OCLH::active) {
		if (opts.getCLDevice() == 0) {
			OCLH::listCLDevices();
			return EXIT_SUCCESS;
		} else {
			try {
				OCLH::prepareCLDevice(opts.getCLDevice()-1);
			} catch (cl::Error& e) {
				Display::printText(e.what());
				Display::printText("Will fall back to sequential version.");
				OCLH::active = false;
			}
		}
	}
	#endif // INOVESA_USE_CL

	size_t nParticles = UINT16_MAX;
	sstream.str("");
	sstream << nParticles;

	Display::printText("Generating initial particle distribution (using "+
					   sstream.str()+" patricles).");
	PhaseSpace* mesh;
	meshindex_t ps_size = 512;
	constexpr double qmax = 5.0;
	constexpr double pmax = 5.0;
	mesh = new PhaseSpace(ps_size,-qmax,qmax,-pmax,pmax);

	std::random_device seed;
	std::default_random_engine engine(seed());

	std::normal_distribution<> x(0.0,1.0);
	std::normal_distribution<> y(0.0,1.0);

	float amplitude = 2.0f;
	float pulselen = 1.90e-3f;
	float wavelen = 6.42e-5f;

	for (unsigned int i=0; i<nParticles; i++) {
		float xf = x(engine);
		float yf = y(engine)
						+ std::exp(-std::pow(xf/(std::sqrt(2)*pulselen/2.35),2))
						* amplitude * std::sin(2*M_PI*xf/wavelen);
		meshindex_t x = std::lround((xf/qmax+0.5f)*ps_size);
		meshindex_t y = std::lround((yf/pmax+0.5f)*ps_size);
		if (x < ps_size && y < ps_size) {
			(*mesh)[x][y] += 2*M_PI*ps_size*ps_size/(qmax*pmax)/nParticles;
		}
	}

	HDF5File file(opts.getOutFile(),ps_size);
	Display::printText("Will save reults to: \""+opts.getOutFile()+'\"');

	PhaseSpace mesh_rotated(*mesh);

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

	Display::printText("Building FokkerPlanckMap.");
	FokkerPlanckMap fpm(&mesh_rotated,mesh,ps_size,ps_size,
						FokkerPlanckMap::FPType::full,e0,
						FokkerPlanckMap::DerivationType::cubic);

	Display::printText("Building RotationMap.");
	RotationMap rm(mesh,&mesh_rotated,ps_size,ps_size,angle,
				   HeritageMap::InterpolationType::cubic,
				   RotationMap::RotationCoordinates::norm_pm1,true);

	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		mesh->syncCLMem(vfps::PhaseSpace::clCopyDirection::cpu2dev);
	}
	#endif // INOVESA_USE_CL
	Display::printText("Starting the simulation.");
	for (unsigned int i=0;i<=steps*rotations;i++) {
		if (i%outstep == 0) {
			#ifdef INOVESA_USE_CL
			if (OCLH::active) {
				mesh->syncCLMem(vfps::PhaseSpace::clCopyDirection::dev2cpu);
			}
			#endif // INOVESA_USE_CL
			file.append(mesh);
			#ifdef INOVESA_USE_GUI
			if (opts.showPhaseSpace()) {
				display->createTexture(mesh);
				display->draw();
				display->delTexture();
			} else
			#endif
			{
				std::stringstream status;
				status << static_cast<float>(i)/steps << '/' << rotations;
				Display::printText(status.str());
			}
		}
		rm.apply();
		fpm.apply();
	}
	#ifdef INOVESA_USE_CL
	if (OCLH::active) {
		OCLH::queue.flush();
	}
	#endif // INOVESA_USE_CL

	delete mesh;

	Display::printText("Finished.");

	return EXIT_SUCCESS;
}

