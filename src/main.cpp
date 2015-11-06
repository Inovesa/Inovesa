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
#include "IO/GUI/Plot2DLine.hpp"
#include "IO/GUI/Plot3DColormap.hpp"
#include "PhaseSpace.hpp"
#include "Impedance.hpp"
#include "CL/OpenCLHandler.hpp"
#include "HM/FokkerPlanckMap.hpp"
#include "HM/Identity.hpp"
#include "HM/KickMap.hpp"
#include "HM/RotationMap.hpp"
#include "HM/WakeFunctionMap.hpp"
#include "HM/WakeImpedanceMap.hpp"
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

    opts.save("foo.cfg");

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
    const double qmin = -qmax;
    const double pmin = qmin;
    const double pmax = qmax;

    std::string startdistfile = opts.getStartDistFile();

    if (startdistfile.length() <= 4) {
        ps_size = opts.getMeshSize();
        if (ps_size == 0) {
            Display::printText("Please give file for initial distribution "
                               "or size of target mesh > 0.");
        }
        Display::printText("Generating (gaussian) initial distribution.");
        mesh = new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,
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

            mesh = new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,
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
        ps_size = opts.getMeshSize();
        mesh = new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,
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

    double bl = opts.getNaturalBunchLength();
    double f0 = opts.getRevolutionFrequency();
    unsigned int padding =opts.getPadding();

    Impedance* impedance = nullptr;
    if (opts.getImpedanceFile() == "") {
        Display::printText("Will use free space CSR impedance. "
                           "(Give impedance file for other impedance model.)");
        impedance = new Impedance(Impedance::ImpedanceModel::FreeSpaceCSR,
                                  ps_size*std::max(padding,1u),f0,
                                  ps_size*vfps::physcons::c/(2*qmax*bl),false);
    } else {
        Display::printText("Reading impedance from: \""
                           +opts.getImpedanceFile()+"\"");
        impedance = new Impedance(opts.getImpedanceFile(),
                                  ps_size*vfps::physcons::c/(2*qmax*bl));
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
    WakeKickMap* wkm = nullptr;
    WakeFunctionMap* wfm = nullptr;
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
        Display::printText("Building WakeFunctionMap.");
        wfm = new WakeFunctionMap(mesh_damdiff,mesh,ps_size,ps_size,
                                  wake,HeritageMap::InterpolationType::cubic);
        wkm = wfm;
    } else {
        double Ib = opts.getBunchCurrent();
        double E0 = opts.getBeamEnergy();
        double sigmaE = opts.getEnergySpread();
        double rb = opts.getBendingRadius();
        Display::printText("Calculating WakeFunction.");
        field = new ElectricField(mesh,impedance,padding*ps_size,padding,true);
        Display::printText("Building WakeFunctionMap.");
        wkm = new WakeImpedanceMap(mesh_damdiff,mesh,ps_size,ps_size,
                                   field,impedance,
                                   HeritageMap::InterpolationType::cubic);
    }

    HDF5File* file = nullptr;
    std::string ofname = opts.getOutFile();
    if ( isOfFileType(".h5",ofname)) {
        std::string cfgname = ofname.substr(0,ofname.find(".h5"))+".cfg";
        opts.save(cfgname);
        Display::printText("Saved configuiration to: \""+cfgname+'\"');
        file = new HDF5File(ofname,mesh,field,impedance,wfm);
        Display::printText("Will save results to: \""+ofname+'\"');
    } else {
        Display::printText("Information: Will not save results.");
    }

    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        mesh->syncCLMem(vfps::PhaseSpace::clCopyDirection::cpu2dev);
    }
    #endif // INOVESA_USE_CL
    #ifdef INOVESA_USE_GUI
    Plot2DLine* bpv = nullptr;
    Plot3DColormap* psv = nullptr;
    bool gui = opts.showPhaseSpace();
    if (gui) {
        try {
            psv = new Plot3DColormap();
            display->addElement(psv);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(psv);
            delete psv;
            psv = nullptr;
            gui = false;
        }
    }
    if (gui) {
        try {
            bpv = new Plot2DLine();
            display->addElement(bpv);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(bpv);
            delete bpv;
            bpv = nullptr;
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
                file->append(wkm);
            }
            #ifdef INOVESA_USE_GUI
            if (gui) {
                psv->createTexture(mesh);
                bpv->createLine(mesh);
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
        file->append(wkm);
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
