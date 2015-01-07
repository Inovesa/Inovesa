#include <fstream>
#include <iostream>
#include <sstream>

#include "Display.hpp"
#include "Mesh2D.hpp"
#include "CL/CLProgs.hpp"
#include "CL/OpenCLHandler.hpp"
#include "IO/HDF5File.hpp"

#include "main.hpp"

enum class rmtype {
	mesh, space, normal
};

enum class pattern {
	square, gaus, half
};

int main(int argc, char** argv)
{
	vfps::Mesh2D<meshdata_t> mesh(vfps::fs_xsize,-10.0,10.0,
								  vfps::fs_ysize,-10.0,10.0);

	vfps::HDF5File file("results.h5");
    
    
    constexpr double rotations = 0.125;
	constexpr unsigned int patterndim_x = 512;
	constexpr unsigned int patterndim_y = 768;

	constexpr rmtype rm = rmtype::space;
    constexpr pattern ptrntype = pattern::half;

	// no more settings below this line

	std::stringstream nrbuf;

	std::string resdir("results/");

	nrbuf.str("");
	nrbuf << vfps::fs_xsize;
	std::string meshdim = nrbuf.str();

	nrbuf.str("");
	nrbuf << patterndim_x;
	std::string pattern_x = nrbuf.str();

	nrbuf.str("");
	nrbuf << patterndim_y;
	std::string pattern_y = nrbuf.str();

	nrbuf.str("");
	nrbuf << rotations*100;
	std::string rotation = nrbuf.str();

	std::string rottype;
	switch (rm) {
	case rmtype::space:
		rottype = "space";
		break;
	case rmtype::normal:
		rottype = "normal";
		break;
	case rmtype::mesh:
	default:
		rottype = "mesh";
		break;
	}

	std::string ptrn;
	switch (ptrntype) {
	case pattern::square:
		ptrn = "sqr";
		break;
	case pattern::gaus:
	default:
		ptrn = "gaus";
		break;
	case pattern::half:
		ptrn = "half";
		break;
	}

	std::string fname = meshdim + "d-" + rottype + "_"
					 + ptrn + "-" + pattern_x+ "-" +pattern_y +"_"
					 + rotation+ ".dat";

	std::ofstream results(resdir+ "final_" + fname);
	std::ofstream decay(resdir+ "decay_" + fname);

	switch (ptrntype) {
	case pattern::square:
		for (unsigned int x=vfps::fs_xsize/4; x<vfps::fs_xsize*3/4; x++) {
			for (unsigned int y=vfps::fs_ysize/4; y<vfps::fs_ysize*3/4; y++) {
				mesh[x][y] = 1.0f;
			}
		}
		for (unsigned int x=vfps::fs_xsize/4; x<vfps::fs_xsize*3/4; x++) {
				mesh[x][vfps::fs_ysize/2] = 0.5f;
		}
		for (unsigned int y=vfps::fs_ysize/4; y<vfps::fs_ysize*3/4; y++) {
				mesh[vfps::fs_xsize/2][y] = 0.5f;
		}
		for (unsigned int x=vfps::fs_xsize/4; x<vfps::fs_xsize*3/4; x++) {
			mesh[x][x] = 0.0f;
			mesh[x][vfps::fs_ysize-x] = 0.0f;
		}
		break;
	case pattern::gaus:
	default:
		for (int x=0; x<int(vfps::fs_xsize); x++) {
			for (int y=0; y<int(vfps::fs_ysize); y++) {
				mesh[x][y] = std::exp(-std::pow(x-int(vfps::fs_xsize)/2,2)/patterndim_x
									  -std::pow(y-int(vfps::fs_ysize)/2,2)/patterndim_y);
			}
		}
		break;
	case pattern::half:
		for (int x=0; x<int(vfps::fs_xsize/2); x++) {
			for (int y=0; y<int(vfps::fs_ysize); y++) {
				mesh[x][y] = 1;
			}
		}
		for (int x=0; x<int(vfps::fs_xsize/2); x++) {
			mesh[x][x] = 0;
			mesh[x][vfps::fs_ysize/2] = 0;
			mesh[x][vfps::fs_ysize-x] = 0;
		}
		break;
	}
	file.write(&mesh);

	unsigned int device = 0;
	if (argc == 2 ) {
		std::stringstream dev(argv[1]);
		dev >> device;
		device--;
	}

#ifdef FR_USE_GUI
	Display display;
    display.createTexture(&mesh);
	display.draw();
#endif

#ifdef FR_USE_CL
	prepareCLEnvironment(device);
	prepareCLProgs();
#else
	std::vector<meshdata_t> kick(vfps::fs_xsize,0.0);
#endif
	constexpr unsigned int steps = 4000;
	constexpr double angle = 2*M_PI/steps;
	switch (rm) {
	case rmtype::space:
        mesh.setRotationMap(angle,vfps::Mesh2D<meshdata_t>::MESH_NORMALIZATION::PHYSICAL);
		break;
	case rmtype::normal:
        mesh.setRotationMap(angle,vfps::Mesh2D<meshdata_t>::MESH_NORMALIZATION::ONE);
		break;
	case rmtype::mesh:
	default:
        mesh.setRotationMap(angle,vfps::Mesh2D<meshdata_t>::MESH_NORMALIZATION::NONE);
		break;
	}
#ifdef FR_USE_CL
	mesh.__initOpenCL();
	mesh.syncData();
#endif
	unsigned int outstep = 1;
	unsigned int i;
	for (i=0;i<steps*rotations;i++) {
        if (i%outstep != 0) {
            mesh.rotate();
		} else {
			if (i == 10)
				outstep = 10;
			if (i == 100)
				outstep = 100;
#ifdef FR_USE_GUI
			display.createTexture(&mesh);
			display.draw();
#endif
			file.append(&mesh);

			meshdata_t sum = 0.0;
			for (unsigned int q_i= floor(vfps::fs_xsize/2.0)-2; q_i < ceil(vfps::fs_xsize/2.0)+2; q_i ++) {
				meshdata_t linesum = 0.0;
				for (unsigned int p_i= floor(vfps::fs_ysize/2.0)-2; p_i < ceil(vfps::fs_ysize/2.0)+2; p_i ++) {
					linesum += mesh[q_i][p_i];
				}
				sum += linesum;
			}

			decay << double(i)/double(steps) << '\t'
				  << mesh[vfps::fs_xsize/2][vfps::fs_ysize/2] - 1.0 << '\t'
				  << sum << std::endl;
			std::cout << double(i)/double(steps) << '\t'
					  << mesh[vfps::fs_xsize/2][vfps::fs_ysize/2] - 1.0 << '\t'
					  << sum << std::endl;
            mesh.rotate();
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
	display.delTexture();
#endif
#ifdef FR_USE_CL
	OCLH::queue.flush();
#endif


	meshdata_t sum = 0.0;
	for (unsigned int q_i= floor(vfps::fs_xsize/2.0)-2; q_i < ceil(vfps::fs_xsize/2.0)+2; q_i ++) {
		meshdata_t linesum = 0.0;
		for (unsigned int p_i= floor(vfps::fs_ysize/2.0)-2; p_i < ceil(vfps::fs_ysize/2.0)+2; p_i ++) {
			linesum += mesh[q_i][p_i];
		}
		sum += linesum;
	}

	decay << double(i)/double(steps) << '\t'
		  << mesh[vfps::fs_xsize/2][vfps::fs_ysize/2] -1.0 << '\t'
		  << sum << std::endl;
	std::cout << double(i)/double(steps) << '\t'
			  << mesh[vfps::fs_xsize/2][vfps::fs_ysize/2] - 1.0 << '\t'
			  << sum << std::endl;

	for (unsigned int x=0; x<vfps::fs_xsize; x++) {
		for (unsigned int y=0; y<vfps::fs_ysize; y++) {
			results << mesh[x][y] << '\t';
		}
		results << std::endl;
    }

	return EXIT_SUCCESS;
}

