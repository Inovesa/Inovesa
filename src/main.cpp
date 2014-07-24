#include <fstream>
#include <iostream>

#include "Display.hpp"
#include "Mesh2D.hpp"

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

	Display display;

	std::vector<meshdata_t> kick(xsize,0.0);
	for (unsigned int i=0;i<1000;i++) {
		mesh.rotate(0.01);
		//mesh.rotateAndKick(0.01,kick);
		display.draw(&mesh);
	}

    return 0;
}

