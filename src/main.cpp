#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Display.hpp"
#include "PhaseSpace.hpp"
#include "CL/CLProgs.hpp"
#include "CL/OpenCLHandler.hpp"
#include "IO/HDF5File.hpp"

#include "main.hpp"

enum class pattern {
	square, gaus, half, quarters
};

int main(int argc, char** argv)
{
	vfps::PhaseSpace mesh(	vfps::ps_xsize,-10.0,10.0,
							vfps::ps_ysize,-10.0,10.0);

	vfps::HDF5File file("results.h5");
    
    
	constexpr double rotations = 1;
	constexpr unsigned int patterndim_x = 512;
	constexpr unsigned int patterndim_y = 4048;
	constexpr unsigned int pattern_margin = 128;
	constexpr unsigned int pattern_max = 1.0f;

	constexpr vfps::PhaseSpace::ROTATION_TYPE rt
			= vfps::PhaseSpace::ROTATION_TYPE::MESH;
	constexpr pattern ptrntype = pattern::quarters;

	// no more settings below this line

	std::stringstream nrbuf;

	std::string resdir("results/");

	nrbuf.str("");
	nrbuf << vfps::ps_xsize;
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
	switch (rt) {
	case vfps::PhaseSpace::ROTATION_TYPE::SPACE:
		rottype = "space";
		break;
	case vfps::PhaseSpace::ROTATION_TYPE::NORMAL:
		rottype = "normal";
		break;
	case vfps::PhaseSpace::ROTATION_TYPE::MESH:
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
		for (unsigned int x=vfps::ps_xsize/4; x<vfps::ps_xsize*3/4; x++) {
			for (unsigned int y=vfps::ps_ysize/4; y<vfps::ps_ysize*3/4; y++) {
				mesh[x][y] = 1.0f*pattern_max;
			}
		}
		for (unsigned int x=vfps::ps_xsize/4; x<vfps::ps_xsize*3/4; x++) {
				mesh[x][vfps::ps_ysize/2] = 0.5f*pattern_max;
		}
		for (unsigned int y=vfps::ps_ysize/4; y<vfps::ps_ysize*3/4; y++) {
				mesh[vfps::ps_xsize/2][y] = 0.5f*pattern_max;
		}
		for (unsigned int x=vfps::ps_xsize/4; x<vfps::ps_xsize*3/4; x++) {
			mesh[x][x] = 0.0f;
			mesh[x][vfps::ps_ysize-x] = 0.0f;
		}
		break;
	case pattern::gaus:
	default:
		for (int x=0; x<int(vfps::ps_xsize); x++) {
			for (int y=0; y<int(vfps::ps_ysize); y++) {
				mesh[x][y] = std::exp(-std::pow(x-int(vfps::ps_xsize)/2,2)/patterndim_x
									  -std::pow(y-int(vfps::ps_ysize)/2,2)/patterndim_y);
			}
		}
		break;
	case pattern::half:
		for (int y=pattern_margin; y<int(vfps::ps_ysize-pattern_margin); y++) {
			for (int x=pattern_margin; x<int(vfps::ps_xsize/2); x++) {
				mesh[x][y] = 1.0f*pattern_max;
			}
		}
		for (int x=pattern_margin; x<int(vfps::ps_xsize/2); x++) {
			mesh[x][x] = 0;
			mesh[x][vfps::ps_ysize/2] = 0;
			mesh[x][vfps::ps_ysize-x] = 0;
		}
		break;
	case pattern::quarters:
		for (int y=pattern_margin; y<int(vfps::ps_ysize/2); y++) {
			for (int x=pattern_margin; x<int(vfps::ps_xsize/2); x++) {
				mesh[x][y] = (float(y-pattern_margin))
							/float(vfps::ps_ysize/2-pattern_margin)*pattern_max;
			}
		}
		for (int y=vfps::ps_ysize/2; y<int(vfps::ps_ysize-pattern_margin); y++) {
			for (int x=pattern_margin; x<int(vfps::ps_xsize/2); x++) {
				mesh[x][y] = 1.0f*pattern_max;
			}
		}
		for (int x=pattern_margin; x<int(vfps::ps_xsize/2); x++) {
			mesh[x][vfps::ps_ysize-x] = 0.0f*pattern_max;
		}
		for (int x=vfps::ps_xsize/2; x<int(vfps::ps_xsize-pattern_margin); x++) {
			for (int y=pattern_margin; y<int(vfps::ps_ysize/2); y++) {
				mesh[x][y] = std::exp(-std::pow(x-int(vfps::ps_xsize)/2,2)/patterndim_x
									  -std::pow(y-int(vfps::ps_ysize)/2,2)/patterndim_y)
						*pattern_max;
			}
		}

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
	std::vector<vfps::meshdata_t> kick(vfps::ps_xsize,0.0);
#endif
	constexpr unsigned int steps = 4000;
	constexpr double angle = 2*M_PI/steps;
	mesh.setRotationMap(angle,rt);
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

			vfps::meshdata_t sum = 0.0;
			for (unsigned int q_i= floor(vfps::ps_xsize/2.0)-2; q_i < ceil(vfps::ps_xsize/2.0)+2; q_i ++) {
				vfps::meshdata_t linesum = 0.0;
				for (unsigned int p_i= floor(vfps::ps_ysize/2.0)-2; p_i < ceil(vfps::ps_ysize/2.0)+2; p_i ++) {
					linesum += mesh[q_i][p_i];
				}
				sum += linesum;
			}

			decay << double(i)/double(steps) << '\t'
				  << mesh[vfps::ps_xsize/2][vfps::ps_ysize/2] - 1.0 << '\t'
				  << sum << std::endl;
			#ifdef FR_PRINT_RESULTS
			std::cout << double(i)/double(steps) << '\t'
					  << mesh[vfps::ps_xsize/2][vfps::ps_ysize/2] - 1.0 << '\t'
					  << sum << std::endl;
			#endif
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
	file.append(&mesh);
	display.delTexture();
#endif
#ifdef FR_USE_CL
	OCLH::queue.flush();
#endif


	vfps::meshdata_t sum = 0.0;
	for (unsigned int q_i= floor(vfps::ps_xsize/2.0)-2; q_i < ceil(vfps::ps_xsize/2.0)+2; q_i ++) {
		vfps::meshdata_t linesum = 0.0;
		for (unsigned int p_i= floor(vfps::ps_ysize/2.0)-2; p_i < ceil(vfps::ps_ysize/2.0)+2; p_i ++) {
			linesum += mesh[q_i][p_i];
		}
		sum += linesum;
	}

	decay << double(i)/double(steps) << '\t'
		  << mesh[vfps::ps_xsize/2][vfps::ps_ysize/2] -1.0 << '\t'
		  << sum << std::endl;
	#ifdef FR_PRINT_RESULTS
	std::cout << double(i)/double(steps) << '\t'
			  << mesh[vfps::ps_xsize/2][vfps::ps_ysize/2] - 1.0 << '\t'
			  << sum << std::endl;
	#endif

	for (unsigned int x=0; x<vfps::ps_xsize; x++) {
		for (unsigned int y=0; y<vfps::ps_ysize; y++) {
			results << mesh[x][y] << '\t';
		}
		results << std::endl;
    }

	return EXIT_SUCCESS;
}

