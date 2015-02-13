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

enum class pattern {
	square, gaus, half, quarters
};

using namespace vfps;

int main(int argc, char** argv)
{
	// settings
	constexpr double rotations = 1;
	constexpr unsigned int patterndim_x = 512;
	constexpr unsigned int patterndim_y = 4048;
	constexpr unsigned int pattern_margin = 128;
	float pattern_max;
	if (std::is_same<meshdata_t,fixp32>::value) {
		pattern_max = 1.0;
	} else {
		pattern_max = 0.5;
	}

	constexpr pattern ptrntype = pattern::quarters;
	constexpr unsigned int steps = 4000;

	/* @todo: remove global settings from main.hpp
	 * (HDF5File::HDF5File() could take mesh as argument)
	 */
	PhaseSpace mesh(-10.0,10.0,-10.0,10.0);
	PhaseSpace mesh_rotated(mesh);

	HDF5File file("results.h5");

	// create pattern to start with
	switch (ptrntype) {
	case pattern::square:
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
		break;
	case pattern::gaus:
	default:
		for (int x=0; x<int(ps_xsize); x++) {
			for (int y=0; y<int(ps_ysize); y++) {
				mesh[x][y] = std::exp(-std::pow(x-int(ps_xsize)/2,2)/patterndim_x
									  -std::pow(y-int(ps_ysize)/2,2)/patterndim_y);
			}
		}
		break;
	case pattern::half:
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
		break;
	case pattern::quarters:
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

	}
	file.write(&mesh);

#ifdef FR_USE_GUI
	Display display;
    display.createTexture(&mesh);
	display.draw();
#endif

#ifdef FR_USE_CL
	// OpenCL device can be given as command line argument
	unsigned int device = 0;
	if (argc == 2 ) {
		std::stringstream dev(argv[1]);
		dev >> device;
		device--;
	}
	prepareCLEnvironment(device);
	prepareCLProgs();
#endif
	// angle of one rotation step (in rad)
	constexpr double angle = 2*M_PI/steps;
	RotationMap rm(&mesh,&mesh_rotated,ps_xsize,ps_ysize,angle);
#ifdef FR_USE_CL
	mesh.__initOpenCL();
	mesh.syncData();
#endif
	unsigned int outstep = 100;
	unsigned int i;
	for (i=0;i<steps*rotations;i++) {
        if (i%outstep != 0) {
			rm.apply();
			swap(mesh,mesh_rotated);
		} else {
#ifdef FR_USE_GUI
			display.createTexture(&mesh);
			display.draw();
#endif
			file.append(&mesh);

			rm.apply();
			swap(mesh,mesh_rotated);
#ifdef FR_USE_GUI
			display.delTexture();
#endif
#ifdef FR_USE_CL
			mesh.syncData();
#endif
		}
	}
#ifdef FR_USE_CL
	mesh.syncData();
#endif
#ifdef FR_USE_GUI
	display.createTexture(&mesh);
	display.draw();
	file.append(&mesh);
	display.delTexture();
#endif
#ifdef FR_USE_CL
	OCLH::queue.flush();
#endif

	return EXIT_SUCCESS;
}

