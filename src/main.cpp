/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2015: Patrik Schönfeldt                                 *
 *                                                                            *
 * This file is part of Inovesa.                                              *
 * Inovesa is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Inovesa is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           *
 ******************************************************************************/

#include <chrono>
#include <climits>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#ifdef INOVESA_USE_PNG
#include <png++/png.hpp>
#endif
#include <sstream>

#include "defines.hpp"
#include "IO/Display.hpp"
#include "IO/GUI/Plot2D.hpp"
#include "PhaseSpace.hpp"
#include "Impedance.hpp"
#include "CL/OpenCLHandler.hpp"
#include "HM/FokkerPlanckMap.hpp"
#include "HM/Identity.hpp"
#include "HM/KickMap.hpp"
#include "HM/RotationMap.hpp"
#include "HM/WakeKickMap.hpp"
#include "IO/HDF5File.hpp"
#include "IO/ProgramOptions.hpp"

using namespace vfps;

inline bool isOfFileType(std::string ending, std::string fname)
{
    return ( fname.size() > ending.size() &&
        std::equal(ending.rbegin(), ending.rend(),fname.rbegin()));
}

int main(int argc, char** argv)
{
    Display::start_time = std::chrono::system_clock::now();

    std::time_t start_ctime
            = std::chrono::system_clock::to_time_t(Display::start_time);
    std::stringstream sstream;
    sstream << std::ctime(&start_ctime);

    ProgramOptions opts;

    try {
        if (!opts.parse(argc,argv)) {
            return EXIT_SUCCESS;
        }
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    std::string timestring = sstream.str();
    timestring.resize(timestring.size()-1);

    sstream.str("");
    sstream << INOVESA_VERSION_RELEASE << '.'
            << INOVESA_VERSION_MINOR << '.'
            << INOVESA_VERSION_FIX;

    Display::printText("Started Inovesa (v"
                       +sstream.str()+") at "+timestring);

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

    PhaseSpace* mesh;
    meshindex_t ps_size;
    const double qmax = opts.getPhaseSpaceSize();
    const double pmax = qmax;

    std::string startdistfile = opts.getStartDistFile();

    if (startdistfile.length() <= 4) {
        ps_size = opts.getMeshSize();
        if (ps_size == 0) {
            Display::printText("Please give file for initial distribution "
                               "or size of target mesh > 0.");
        }
        Display::printText("Generating (gaussian) initial distribution.");
        mesh = new PhaseSpace(ps_size,-qmax,qmax,-pmax,pmax,
                              opts.getNaturalBunchLength());
        for (meshindex_t x = 0; x < ps_size; x++) {
            for (meshindex_t y = 0; y < ps_size; y++) {
                (*mesh)[x][y]
                    = std::exp(-std::pow((float(x)/ps_size-0.5f)*2*qmax,2.0f)
                               /2.0f)
                    * std::exp(-std::pow((float(y)/ps_size-0.5f)*2*pmax,2.0f)
                               /2.0f);
            }
        }
    } else {
        Display::printText("Reading in initial distribution from: \""
                           +startdistfile+'\"');
    #ifdef INOVESA_USE_PNG
    // check for file ending .png
    if (isOfFileType(".png",startdistfile)) {
        // load pattern to start with
        png::image<png::gray_pixel_16> image;
        try {
            image.read(opts.getStartDistFile());
        } catch ( const png::std_error &e ) {
            std::cerr << e.what() << std::endl;
            return EXIT_SUCCESS;
        }
        catch ( const png::error &e ) {
            std::cerr << "Problem loading " << startdistfile
                      << ": " << e.what() << std::endl;
            return EXIT_SUCCESS;
        }

        if (image.get_width() == image.get_height()) {
            ps_size = image.get_width();

            mesh = new PhaseSpace(ps_size,-qmax,qmax,-pmax,pmax,
                                  opts.getNaturalBunchLength());

            for (unsigned int x=0; x<ps_size; x++) {
                for (unsigned int y=0; y<ps_size; y++) {
                    (*mesh)[x][y] = image[ps_size-y-1][x]/float(UINT16_MAX);
                }
            }
            std::stringstream imgsize;
            imgsize << ps_size;
            Display::printText("Read phase space (a="+imgsize.str()+" px).");
        } else {
            std::cerr << "Phase space has to be quadratic. Please adjust "
                      << startdistfile << std::endl;

            return EXIT_SUCCESS;
        }
    } else
    #endif // INOVESA_USE_PNG
    if (isOfFileType(".txt",startdistfile)) {
        ps_size = 512;
        mesh = new PhaseSpace(ps_size,-qmax,qmax,-pmax,pmax,
                              opts.getNaturalBunchLength());

        std::ifstream ifs;
        ifs.open(startdistfile);

        ifs.unsetf(std::ios_base::skipws);

        // count the newlines with an algorithm specialized for counting:
        size_t line_count = std::count(
            std::istream_iterator<char>(ifs),
            std::istream_iterator<char>(),
            '\n');

        ifs.setf(std::ios_base::skipws);
        ifs.clear();
        ifs.seekg(0,ifs.beg);

        while (ifs.good()) {
            float xf,yf;
            ifs >> xf >> yf;
            meshindex_t x = std::lround((xf/qmax+0.5f)*ps_size);
            meshindex_t y = std::lround((yf/pmax+0.5f)*ps_size);
            if (x < ps_size && y < ps_size) {
                (*mesh)[x][y] += 1.0/line_count;
            }
        }
        ifs.close();

        // normalize to higest peak
        meshdata_t maxval = std::numeric_limits<meshdata_t>::min();
        for (unsigned int x=0; x<ps_size; x++) {
            for (unsigned int y=0; y<ps_size; y++) {
                if ((*mesh)[x][y] > maxval) {
                    maxval = (*mesh)[x][y];
                }
            }
        }
        for (unsigned int x=0; x<ps_size; x++) {
            for (unsigned int y=0; y<ps_size; y++) {
                (*mesh)[x][y] /= maxval;
            }
        }
    } else {
        Display::printText("Unknown format of input file. Will now quit.");
        return EXIT_SUCCESS;
    }
    }

    Impedance* impedance = nullptr;
    if (opts.getImpedanceFile() == "") {
        Display::printText("No impedance file given. "
                           "Will use free space impedance.");
        impedance = new Impedance(Impedance::ImpedanceModel::FreeSpace,
                                  ps_size*std::max(opts.getPadding(),1u));
    } else {
        Display::printText("Reading impedance from: \""
                           +opts.getImpedanceFile()+"\"");
        impedance = new Impedance(opts.getImpedanceFile());
        if (impedance->maxN() < ps_size) {
            Display::printText("No valid impedance file. "
                               "Will now quit.");
            return EXIT_SUCCESS;
        }
    }

    PhaseSpace* mesh_rotated = new PhaseSpace(*mesh);
    PhaseSpace* mesh_damdiff = new PhaseSpace(*mesh);

    #ifdef INOVESA_USE_GUI
    Display* display = nullptr;
    if (opts.showPhaseSpace()) {
        display = new Display();
    }
    #endif

    const unsigned int steps = std::max(opts.getSteps(),1u);
    const unsigned int outstep = opts.getOutSteps();
    const float rotations = opts.getNRotations();
    const double f_s = opts.getSyncFreq();
    const double t_d = opts.getDampingTime();
    const double dt = 1.0/(f_s*steps);

    HeritageMap* rm;
    if (steps > 1) {
        /* angle of one rotation step (in rad)
         * (angle = 2*pi corresponds to 1 synchrotron period)
         */
        const double angle = 2*M_PI/steps;

        size_t rotmapsize;
        #ifdef INOVESA_USE_CL
        if (OCLH::active) {
            rotmapsize = 0;
        } else
        #endif
        {
            rotmapsize = ps_size*ps_size/2;
        }
        Display::printText("Building RotationMap.");
        rm = new RotationMap(mesh,mesh_rotated,ps_size,ps_size,angle,
                             HeritageMap::InterpolationType::cubic,
                             RotationMap::RotationCoordinates::norm_pm1,
                             true,rotmapsize);
    } else {
        rm = new Identity(mesh,mesh_rotated,ps_size,ps_size);
    }

    double e0;
    if (t_d > 0) {
        e0 = 2.0/(f_s*t_d*steps);
    } else {
        e0=0;
    }

    HeritageMap* fpm;
    if (e0 > 0) {
        Display::printText("Building FokkerPlanckMap.");
        fpm = new FokkerPlanckMap( mesh_rotated,mesh_damdiff,ps_size,ps_size,
                                   FokkerPlanckMap::FPType::full,e0,
                                   FokkerPlanckMap::DerivationType::cubic);
    } else {
        fpm = new Identity(mesh_rotated,mesh_damdiff,ps_size,ps_size);
    }

    ElectricField* field = nullptr;
    HeritageMap* wkm = nullptr;
    std::vector<std::pair<meshaxis_t,double>> wake;
    std::string wakefile = opts.getWakeFile();
    if (wakefile.size() > 4) {
        field = new ElectricField(mesh,impedance);
        Display::printText("Reading WakeFunction from "+wakefile+".");
        std::ifstream ifs;
        ifs.open(wakefile);

        while (ifs.good()) {
            double q,f;
            ifs >> q >> f;
            wake.push_back(std::pair<meshaxis_t,double>(q,f));
        }
        ifs.close();
        Display::printText("Building WakeKickMap.");
        wkm = new WakeKickMap(mesh_damdiff,mesh,ps_size,ps_size,
                              wake,WakeKickMap::InterpolationType::cubic);
    } else {
        double Ib = opts.getBunchCurrent();
        double bl = opts.getNaturalBunchLength();
        double E0 = 1.3e9;
        double sigmaE = 4.7e-4;
        double f0 = opts.getRevolutionFrequency();
        double rb = opts.getBendingRadius();
        Display::printText("Calculating WakeFunction.");
        field = new ElectricField(mesh,impedance,Ib,bl,E0,sigmaE,f_s,f0,dt,rb,2048);
        Display::printText("Building WakeKickMap.");
        wkm = new WakeKickMap(mesh_damdiff,mesh,ps_size,ps_size,
                              field,WakeKickMap::InterpolationType::cubic);
    }

    HDF5File* file = nullptr;
    if ( isOfFileType(".h5",opts.getOutFile())) {
        file = new HDF5File(opts.getOutFile(),mesh,field,impedance,
                            static_cast<WakeKickMap*>(wkm));
        Display::printText("Will save results to: \""+opts.getOutFile()+'\"');
    } else {
        Display::printText("Will not save results.");
    }

    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        mesh->syncCLMem(vfps::PhaseSpace::clCopyDirection::cpu2dev);
    }
    #endif // INOVESA_USE_CL
    #ifdef INOVESA_USE_GUI
    Plot2D* psv = nullptr;
    if (opts.showPhaseSpace()) {
        try {
            psv = new Plot2D();
            display->addElement(psv);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(psv);
            delete psv;
            psv = nullptr;
        }
    }
    #endif

    Display::printText("Starting the simulation.");
    for (unsigned int i=0;i<steps*rotations;i++) {
        if (i%outstep == 0) {
            #ifdef INOVESA_USE_CL
            if (OCLH::active) {
                mesh->syncCLMem(vfps::PhaseSpace::clCopyDirection::dev2cpu);
            }
            #endif // INOVESA_USE_CL
            if (file != nullptr) {
                mesh->integral();
                file->append(mesh);
                field->updateCSRSpectrum();
                file->append(field);
            }
            #ifdef INOVESA_USE_GUI
            if (psv != nullptr) {
                psv->createTexture(mesh);
                display->draw();
                psv->delTexture();
            }
            #endif
            std::stringstream status;
            status << static_cast<float>(i)/steps << '/' << rotations;
            Display::printText(status.str(),2.0f);
        }
        rm->apply();
        fpm->apply();
        wkm->apply();
    }

    // save final result
    if (file != nullptr) {
        mesh->integral();
        file->append(mesh);
        field->updateCSRSpectrum();
        file->append(field);
    }
    if (!opts.showPhaseSpace()) {
        std::stringstream status;
        status << rotations << '/' << rotations;
        Display::printText(status.str());
    }

    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        OCLH::queue.flush();
    }
    #endif // INOVESA_USE_CL

    #ifdef INOVESA_USE_GUI
    delete display;
    delete psv;
    #endif

    delete mesh;
    delete mesh_rotated;
    delete mesh_damdiff;

    delete field;
    delete impedance;

    delete rm;
    delete wkm;
    delete fpm;

    Display::printText("Finished.");

    return EXIT_SUCCESS;
}
