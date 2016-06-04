/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Algorithms   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
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
#include <cmath>
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
#include "impedances/FreeSpaceCSR.hpp"
#include "impedances/ParallelPlatesCSR.hpp"
#include "CL/OpenCLHandler.hpp"
#include "SM/FokkerPlanckMap.hpp"
#include "SM/Identity.hpp"
#include "SM/KickMap.hpp"
#include "SM/DriftMap.hpp"
#include "SM/RotationMap.hpp"
#include "SM/RFKickMap.hpp"
#include "SM/WakeFunctionMap.hpp"
#include "SM/WakePotentialMap.hpp"
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
    sstream << 'v' << INOVESA_VERSION_RELEASE << '.'
            << INOVESA_VERSION_MINOR << '.'
            << INOVESA_VERSION_FIX;
    if (std::string(GIT_BRANCH) != "master") {
        sstream << ", Branch: "<< GIT_BRANCH;
    }

    #ifdef INOVESA_USE_CL
    if (opts.getCLDevice() >= 0)
    #endif // INOVESA_USE_CL
    {
    Display::printText("Started Inovesa ("
                       +sstream.str()+") at "+timestring);
    }



    #ifdef INOVESA_USE_GUI
    bool gui = opts.showPhaseSpace();
    Display* display = nullptr;
    if (gui && opts.getCLDevice() >= 0) {
        try {
        display = new Display(opts.getOpenGLVersion());
        } catch (std::exception& e) {
            std::string msg("ERROR: ");
            Display::printText(msg+e.what());
            delete display;
            display = nullptr;
            gui = false;
        }
    }
    #endif

    #ifdef INOVESA_USE_CL
    OCLH::active = (opts.getCLDevice() != 0);
    if (OCLH::active) {
        try {
            OCLH::prepareCLEnvironment(gui);
        } catch (cl::Error& e) {
            Display::printText(e.what());
            Display::printText("Will fall back to sequential version.");
            OCLH::active = false;
        }
    }
    if (OCLH::active) {
        if (opts.getCLDevice() < 0) {
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

    vfps::FokkerPlanckMap::DerivationType derivationtype;
    switch (opts.getDerivationType()) {
    case 3:
        derivationtype = vfps::FokkerPlanckMap::DerivationType::two_sided;
        break;
    case 4:
    default:
        derivationtype = vfps::FokkerPlanckMap::DerivationType::cubic;
        break;
    }

    vfps::SourceMap::InterpolationType interpolationtype;
    switch (opts.getInterpolationPoints()) {
    case 1:
        interpolationtype = vfps::SourceMap::InterpolationType::none;
        break;
    case 2:
        interpolationtype = vfps::SourceMap::InterpolationType::linear;
        break;
    case 3:
        interpolationtype = vfps::SourceMap::InterpolationType::quadratic;
        break;
    default:
    case 4:
        interpolationtype = vfps::SourceMap::InterpolationType::cubic;
        break;
    }

    bool interpol_clamp = opts.getInterpolationBound();
    bool verbose = opts.getVerbosity();

    PhaseSpace* mesh1;
    meshindex_t ps_size = opts.getMeshSize();
    const double pqsize = opts.getPhaseSpaceSize();
    const double qcenter = -opts.getPSShiftX()*pqsize/(ps_size-1);
    const double pcenter = -opts.getPSShiftY()*pqsize/(ps_size-1);
    const double pqhalf = pqsize/2;
    const double qmax = qcenter + pqhalf;
    const double qmin = qcenter - pqhalf;
    const double pmax = pcenter + pqhalf;
    const double pmin = pcenter - pqhalf;

    const double sE = opts.getEnergySpread();
    const double E0 = opts.getBeamEnergy();
    const double dE = sE*E0;
    const double f0 = opts.getRevolutionFrequency();
    const double R_tmp = opts.getBendingRadius();
    if (R_tmp > 0) {
        Display::printText("Non iso-magneting rings to be implemented.");
    }
    const double R = (R_tmp>0) ? R_tmp : physcons::c/(2*M_PI*f0);
    #ifdef INOVESA_USE_HDF5
    const double fc = opts.getCutoffFrequency();
    #endif // INOVESA_USE_HDF5
    const double H = opts.getHarmonicNumber();
    const double height = opts.getVacuumChamberHeight();
    const double V = opts.getRFVoltage();
    const double fs = opts.getSyncFreq();
    const double bl = physcons::c*dE/H/std::pow(f0,2.0)/V*fs;
    const double Ib = opts.getBunchCurrent();
    const double Fk = opts.getStartDistParam();

    const unsigned int steps = std::max(opts.getSteps(),1u);
    const unsigned int outstep = opts.getOutSteps();
    const float rotations = opts.getNRotations();
    const double t_d = opts.getDampingTime();
    const double dt = 1.0/(fs*steps);
    const double t_sync = 1.0/fs;

    /* angle of one rotation step (in rad)
     * (angle = 2*pi corresponds to 1 synchrotron period)
     */
    const double angle = 2*M_PI/steps;

    std::string startdistfile = opts.getStartDistFile();

    if (height>=0) {
        double shield;
        if (height==0) {
            shield = 0;
        } else {
            shield = bl*std::sqrt(R)*std::pow(height,-3./2.);
        }

        const double Inorm = physcons::IAlfven/physcons::me*2*M_PI
                           * std::pow(dE*fs/f0,2)/V/H* std::pow(bl/R,1./3.);

        const double Ith = Inorm * (0.5+0.34*shield);

        const double S_csr = Ib/Inorm;

        if (verbose) {
            sstream.str("");
            sstream << std::scientific << shield;
            Display::printText("Shielding parameter: "
                               +sstream.str());
            sstream.str("");
            sstream << std::scientific << S_csr;
            Display::printText("CSR strength: "
                               +sstream.str());
            sstream.str("");
            sstream << std::scientific << Ith;
            Display::printText("BBT-Threshold-Current expected at "
                               +sstream.str()+" A.");
        }
    }

    if (verbose) {
        sstream.str("");
        sstream << std::fixed << 1/dt/f0;
        Display::printText("Doing " +sstream.str()+
                           " simulation steps per revolution period.");

        sstream.str("");
        double rotationoffset = std::tan(angle)*ps_size/2;
        sstream << std::fixed << rotationoffset;
        Display::printText("Maximum rotation offset is "
                           +sstream.str()+" (should be < 1).");
    }

    if (startdistfile.length() <= 4) {
        if (ps_size == 0) {
            Display::printText("Please give file for initial distribution "
                               "or size of target mesh > 0.");
        }
        sstream.str("");
        sstream << Fk;
        Display::printText("Generating initial distribution with F(k)="
                           +sstream.str()+".");
        mesh1 = new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,bl,dE,Fk);
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

            mesh1 = new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,bl);

            for (unsigned int x=0; x<ps_size; x++) {
                for (unsigned int y=0; y<ps_size; y++) {
                    (*mesh1)[x][y] = image[ps_size-y-1][x]/float(UINT16_MAX);
                }
            }
            mesh1->syncCLMem(clCopyDirection::cpu2dev);
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
    #ifdef INOVESA_USE_HDF5
    if (isOfFileType(".h5",startdistfile)) {
        mesh1 = new PhaseSpace(HDF5File::readPhaseSpace(startdistfile,qmax,bl,dE));
        mesh1->syncCLMem(clCopyDirection::cpu2dev);
    } else
    #endif
    if (isOfFileType(".txt",startdistfile)) {
        ps_size = opts.getMeshSize();
        mesh1 = new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,bl);

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
                (*mesh1)[x][y] += 1.0/line_count;
            }
        }
        ifs.close();
        mesh1->syncCLMem(clCopyDirection::cpu2dev);
    } else {
        Display::printText("Unknown format of input file. Will now quit.");
        return EXIT_SUCCESS;
    }
    }

    // normalize integral to 1
    mesh1->normalize();

    // find highest peak (for display)
    meshdata_t maxval = std::numeric_limits<meshdata_t>::min();
    for (unsigned int x=0; x<ps_size; x++) {
        for (unsigned int y=0; y<ps_size; y++) {
            maxval = std::max(maxval,(*mesh1)[x][y]);
        }
    }

    if (verbose) {
    sstream.str("");
    sstream << std::fixed << maxval*Ib/f0/physcons::e;
    Display::printText("Maximum particles per grid cell is "
                       +sstream.str()+".");
    }

    const unsigned int padding =std::max(opts.getPadding(),1u);

    Impedance* impedance = nullptr;
    if (opts.getImpedanceFile() == "") {
        if (height>0) {
            Display::printText("Will use parallel plates CSR impedance.");
            impedance = new ParallelPlatesCSR(ps_size*padding,f0,
                                     ps_size*vfps::physcons::c/(2*qmax*bl),height);
        } else {
            Display::printText("Will use free space CSR impedance.");
            impedance = new FreeSpaceCSR(ps_size*padding,f0,
                                         ps_size*vfps::physcons::c/(2*qmax*bl));
        }
    } else {
        Display::printText("Reading impedance from: \""
                           +opts.getImpedanceFile()+"\"");
        impedance = new Impedance(opts.getImpedanceFile(),
                                  ps_size*vfps::physcons::c/(2*qmax*bl));
        if (impedance->nFreqs() < ps_size) {
            Display::printText("No valid impedance file. "
                               "Will now quit.");
            return EXIT_SUCCESS;
        }
    }

    PhaseSpace* mesh2 = new PhaseSpace(*mesh1);
    PhaseSpace* mesh3 = new PhaseSpace(*mesh1);

    SourceMap* rm1;
    SourceMap* rm2 = nullptr;
    const uint_fast8_t rotationtype = opts.getRotationType();
    switch (rotationtype) {
    case 0:
        Display::printText("Initializing RotationMap.");
        rm1 = new RotationMap(mesh1,mesh3,ps_size,ps_size,angle,
                             interpolationtype,interpol_clamp,
                             RotationMap::RotationCoordinates::norm_pm1,0);
        break;
    case 1:
        Display::printText("Building RotationMap.");
        rm1 = new RotationMap(mesh1,mesh3,ps_size,ps_size,angle,
                             interpolationtype,interpol_clamp,
                             RotationMap::RotationCoordinates::norm_pm1,
                             ps_size*ps_size);
        break;
    case 2:
    default:
        Display::printText("Building RFKickMap.");
        rm1 = new RFKickMap(mesh1,mesh2,ps_size,ps_size,angle,
                            interpolationtype,interpol_clamp);

        Display::printText("Building DriftMap.");
        rm2 = new DriftMap(mesh2,mesh3,ps_size,ps_size,angle,
                           interpolationtype,interpol_clamp);
        break;
    }

    double e1;
    if (t_d > 0) {
        e1 = 1.0/(fs*t_d*steps);
    } else {
        e1=0;
    }

    SourceMap* fpm;
    if (e1 > 0) {
        Display::printText("Building FokkerPlanckMap.");
        fpm = new FokkerPlanckMap( mesh3,mesh2,ps_size,ps_size,
                                   FokkerPlanckMap::FPType::full,e1,
                                   derivationtype);
    } else {
        fpm = new Identity(mesh3,mesh2,ps_size,ps_size);
    }

    ElectricField* field = nullptr;
    SourceMap* wm = nullptr;
    WakeKickMap* wkm = nullptr;
    WakeFunctionMap* wfm = nullptr;
    std::vector<std::pair<meshaxis_t,double>> wake;
    std::string wakefile = opts.getWakeFile();
    if (wakefile.size() > 4) {
        field = new ElectricField(mesh2,impedance);
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
        wfm = new WakeFunctionMap(mesh2,mesh1,ps_size,ps_size,
                                  wake,interpolationtype,interpol_clamp);
        wkm = wfm;
    } else {
        Display::printText("Calculating WakePotential.");
        field = new ElectricField(mesh2,impedance,Ib,E0,sE,dt);
        if (height >= 0) {
            Display::printText("Building WakeKickMap.");
            wkm = new WakePotentialMap(mesh2,mesh1,ps_size,ps_size,field,
                                       interpolationtype,interpol_clamp);
        }
    }
    if (wkm != nullptr) {
        wm = wkm;
    } else {
        wm = new Identity(mesh2,mesh1,ps_size,ps_size);
    }

    std::string ofname = opts.getOutFile();
    #ifdef INOVESA_USE_HDF5
    HDF5File* hdf_file = nullptr;
    if ( isOfFileType(".h5",ofname)) {
        std::string cfgname = ofname.substr(0,ofname.find(".h5"))+".cfg";
        opts.save(cfgname);
        Display::printText("Saved configuiration to \""+cfgname+"\".");
        hdf_file = new HDF5File(ofname,mesh1,field,impedance,wfm,Ib,t_sync);
        Display::printText("Will save results to \""+ofname+"\".");
    } else
    #endif // INOVESA_USE_HDF5
    #ifdef INOVESA_USE_PNG
    if ( isOfFileType(".png",ofname)) {
        std::string cfgname = ofname.substr(0,ofname.find(".png"))+".cfg";
        opts.save(cfgname);
        Display::printText("Saved configuiration to \""+cfgname+"\".");
        Display::printText("Will save results to \""+ofname+"\".");
    } else
    #endif // INOVESA_USE_PNG
    {
        Display::printText("Will not save results.");
    }

    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        mesh1->syncCLMem(clCopyDirection::cpu2dev);
    }
    #endif // INOVESA_USE_CL
    #ifdef INOVESA_USE_GUI
    Plot2DLine* bpv = nullptr;
    Plot3DColormap* psv = nullptr;
    Plot2DLine* wpv = nullptr;
    if (gui) {
        try {
            psv = new Plot3DColormap(maxval);
            display->addElement(psv);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(psv);
            delete psv;
            psv = nullptr;
            gui = false;
        }
        try {
            bpv = new Plot2DLine(std::array<float,3>{{1,0,0}});
            display->addElement(bpv);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(bpv);
            delete bpv;
            bpv = nullptr;
        }
        if (wkm != nullptr) {
            try {
                wpv = new Plot2DLine(std::array<float,3>{{0,0,1}});
                display->addElement(wpv);
            } catch (std::exception &e) {
                std::cerr << e.what() << std::endl;
                display->takeElement(wpv);
                delete wpv;
                wpv = nullptr;
            }
        }
    }
    #endif // INOVESSA_USE_GUI

    Display::printText("Starting the simulation.");
    // first update to have wkm for step0 in file
    if (wkm != nullptr) {
        wkm->update();
    }

    for (unsigned int i=0;i<steps*rotations;i++) {
        if (outstep > 0 && i%outstep == 0) {
            #ifdef INOVESA_USE_CL
            if (OCLH::active) {
                mesh1->syncCLMem(clCopyDirection::dev2cpu);
                if (wkm != nullptr) {
                    wkm->syncCLMem(clCopyDirection::dev2cpu);
                }
            }
            #endif // INOVESA_USE_CL
            mesh1->normalize();
            #ifdef INOVESA_USE_HDF5
            if (hdf_file != nullptr) {
                hdf_file->appendTime(static_cast<double>(i)
                                /static_cast<double>(steps));
                mesh1->variance(0);
                mesh1->variance(1);
                hdf_file->append(mesh1);
                field->updateCSR(fc);
                hdf_file->append(field);
                if (wkm != nullptr) {
                    hdf_file->append(wkm);
                }
            }
            #endif // INOVESA_USE_HDF5
            #ifdef INOVESA_USE_GUI
            if (gui) {
                if (psv != nullptr) {
                    psv->createTexture(mesh1);
                }
                if (bpv != nullptr) {
                    bpv->updateLine(mesh1->nMeshCells(0),
                                    mesh1->getProjection(0));
                }
                if (wpv != nullptr) {
                    wpv->updateLine(mesh1->nMeshCells(0),
                                    wkm->getForce());
                }
                display->draw();
                if (psv != nullptr) {
                    psv->delTexture();
                }
            }
            #endif // INOVESSA_USE_GUI
            std::stringstream status;
            status.precision(3);
            status << std::setw(4) << static_cast<float>(i)/steps
                   << '/' << rotations;
            status << "\t1-Q/Q_0=" << 1.0 - mesh1->getIntegral();
            Display::printText(status.str(),2.0f);
        }
        rm1->apply();
        if (rm2 != nullptr) {
          rm2->apply();
        }
        fpm->apply();
        wm->apply();
    }

    #ifdef INOVESA_USE_HDF5
    // save final result
    if (hdf_file != nullptr) {
        hdf_file->appendTime(rotations);
        mesh1->integral();
        hdf_file->append(mesh1);
        field->updateCSR(fc);
        hdf_file->append(field);
        if (wkm != nullptr) {
            hdf_file->append(wkm);
        }
    }
    #endif // INOVESA_USE_HDF5
    #ifdef INOVESA_USE_PNG
    if ( isOfFileType(".png",ofname)) {
        meshdata_t maxval = std::numeric_limits<meshdata_t>::min();
        meshdata_t* val = mesh1->getData();
        for (meshindex_t i=0; i<mesh1->nMeshCells(); i++) {
            maxval = std::max(val[i],maxval);
        }
        png::image< png::gray_pixel_16 > png_file(ps_size, ps_size);
        for (unsigned int x=0; x<ps_size; x++) {
            for (unsigned int y=0; y<ps_size; y++) {
                png_file[ps_size-y-1][x]=(*mesh1)[x][y]/maxval*float(UINT16_MAX);
            }
        }
        png_file.write(ofname);
    }
    #endif

    std::stringstream status;
    status.precision(3);
    status << std::setw(4) << rotations << '/' << rotations;
    status << "\t1-Q/Q_0=" << 1.0 - mesh1->integral();
    Display::printText(status.str());

    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        OCLH::queue.flush();
    }
    #endif // INOVESA_USE_CL

    #ifdef INOVESA_USE_GUI
    delete display;
    delete psv;
    #endif

    delete mesh1;
    delete mesh2;
    delete mesh3;

    delete field;
    delete impedance;

    delete rm1;
    delete rm2;
    delete wm;
    delete fpm;

    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        OCLH::teardownCLEnvironment();
    }
    #endif // INOVESA_USE_CL

    Display::printText("Finished.");

    return EXIT_SUCCESS;
}

