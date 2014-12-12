#include <fstream>
#include <iostream>
#include <sstream>

#include "Display.hpp"
#include "Mesh2D.hpp"
#include "CL/CLProgs.hpp"
#include "CL/OpenCLHandler.hpp"

int main(int argc, char** argv)
{
	constexpr unsigned int xsize = 512;
	constexpr unsigned int ysize = 512;
	vfps::Mesh2D<meshdata_t> mesh(xsize,-10.0,10.0,ysize,-10.0,10.0);

	for (unsigned int x=0; x<xsize/2; x++) {
		for (unsigned int y=0; y<ysize; y++) {
			mesh[x][y] = 1.0f;
		}
	}
	for (unsigned int x=0; x<xsize/2; x++) {
			mesh[x][ysize/2] = 0.5f;
	}
	for (unsigned int x=0; x<xsize/2; x++) {
		mesh[x][x] = 0.0f;
		mesh[x][ysize-x] = 0.0f;
	}

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
	std::vector<meshdata_t> kick(xsize,0.0);
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

