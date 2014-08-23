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
			mesh[x][y] = 1.0;
		}
	}

	unsigned int device = 0;
	if (argc == 2 ) {
		std::stringstream dev(argv[1]);
		dev >> device;
		device--;
	}

	prepareCLEnvironment(device);
	prepareCLProgs();
	Display display;

	std::vector<meshdata_t> kick(xsize,0.0);
	mesh.setRotationMap(0.01);
	mesh.__initOpenCL();
	for (unsigned int i=0;i<1000;i++) {
		display.draw(&mesh);
		mesh.rotate();
		//mesh.rotateAndKick(0.01,kick);
	}

	return EXIT_SUCCESS;
}

