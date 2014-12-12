#include <fstream>
#include <iostream>
#include <sstream>

#include "Display.hpp"
#include "Mesh2D.hpp"
#include "CL/CLProgs.hpp"
#include "CL/OpenCLHandler.hpp"
#include "IO/HDF5File.hpp"

#include "main.hpp"

int main(int argc, char** argv)
{
	vfps::Mesh2D<meshdata_t> mesh(vfps::fs_xsize,-10.0,10.0,
								  vfps::fs_ysize,-10.0,10.0);

	vfps::HDF5File file("results.h5");

	for (unsigned int x=0; x<vfps::fs_xsize/2; x++) {
		for (unsigned int y=0; y<vfps::fs_ysize; y++) {
			mesh[x][y] = 1.0f;
		}
	}
	for (unsigned int x=0; x<vfps::fs_xsize/2; x++) {
			mesh[x][vfps::fs_ysize/2] = 0.5f;
	}
	for (unsigned int x=0; x<vfps::fs_xsize/2; x++) {
		mesh[x][x] = 0.0f;
		mesh[x][vfps::fs_ysize-x] = 0.0f;
	}
	file.write(&mesh);

	unsigned int device = 0;
	if (argc == 2 ) {
		std::stringstream dev(argv[1]);
		dev >> device;
		device--;
	}

	Display display;
	display.createTexture(&mesh);
	display.draw();

#ifdef FR_USE_CL
	prepareCLEnvironment(device);
	prepareCLProgs();
#else
	std::vector<meshdata_t> kick(vfps::fs_xsize,0.0);
#endif
	constexpr unsigned int steps = 1000;
	constexpr float angle = M_PI/2/steps;
#ifdef FR_USE_CL
	mesh.setRotationMap(angle);
	mesh.__initOpenCL(display.getTexture());
	mesh.syncData();
#endif
	for (unsigned int i=0;i<steps;i++) {
		if (i%100 != 0) {
#ifdef FR_USE_CL
			mesh.rotate();
#else
			mesh.rotateAndKick(angle,kick);
#endif
		} else {
			file.append(&mesh);
			display.createTexture(&mesh);
			display.draw();
#ifdef FR_USE_CL
			mesh.rotate();
#else
			mesh.rotateAndKick(angle,kick);
#endif
			display.delTexture();
#ifdef FR_USE_CL
			mesh.syncData();
#endif
		}
	}
	//display.createTexture(&mesh);
	display.draw();
	display.delTexture();
#ifdef FR_USE_CL
	OCLH::queue.flush();
#endif

	return EXIT_SUCCESS;
}

