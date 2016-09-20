/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
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
#include "PS/PhaseSpace.hpp"
#include "Z/FreeSpaceCSR.hpp"
#include "Z/ParallelPlatesCSR.hpp"
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

    std::string ofname = opts.getOutFile();
    #ifdef INOVESA_USE_CL
    if (opts.getCLDevice() >= 0)
    #endif // INOVESA_USE_CL
    {
    if (ofname != "/dev/null") {
        Display::logfile.open(ofname+".log");
    }
    Display::printText("Started Inovesa ("
                       +vfps::inovesa_version()+") at "+timestring);
    if (ofname != "/dev/null") {
        Display::printText("Will create log at \""+ofname+".log\".");
    }
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
    if (opts.getCLDevice() < 0) {
        OCLH::listCLDevices();
        return EXIT_SUCCESS;
    }

    OCLH::active = (opts.getCLDevice() > 0);
    if (OCLH::active) {
        try {
            OCLH::prepareCLEnvironment(gui,opts.getCLDevice()-1);
            std::atexit(OCLH::teardownCLEnvironment);
        } catch (cl::Error& e) {
            Display::printText(e.what());
            Display::printText("Will fall back to sequential version.");
            OCLH::active = false;
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

    const bool interpol_clamp = opts.getInterpolationBound();
    const bool verbose = opts.getVerbosity();
    const bool renormalize = opts.getRenormalizeCharge();

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

    const double sE = opts.getEnergySpread(); // relative energy spread
    const double E0 = opts.getBeamEnergy(); // relative energy spread
    const double dE = sE*E0; // absolute energy spread
    const double f_rev = opts.getRevolutionFrequency();
    const double R_tmp = opts.getBendingRadius();
    const double R_bend = (R_tmp>0) ? R_tmp : physcons::c/(2*M_PI*f_rev);
    const double f0 = (R_tmp<=0) ? f_rev : physcons::c/(2*M_PI*R_bend);

    // scaling for isomagnetic approximation, defined to be <= 1
    const double isoscale = f_rev/f0;

    const double fc = opts.getCutoffFrequency();
    const double H_unscaled = opts.getHarmonicNumber();
    const double H = isoscale*H_unscaled;
    const double gap = opts.getVacuumChamberGap();
    const double V = opts.getRFVoltage();
    const meshaxis_t alpha0 = opts.getAlpha0();
    const meshaxis_t alpha1 = opts.getAlpha1();
    const meshaxis_t alpha2 = opts.getAlpha2();

    // real synchrotron frequency
    const double fs_unscaled = f_rev*std::sqrt(alpha0*H_unscaled*V/(2*M_PI*E0));

    // synchrotron frequency (isomagnetic ring)
    const double fs = fs_unscaled/isoscale;

    // natural RMS bunch length
    const double bl = physcons::c*dE/H/std::pow(f0,2.0)/V*fs;
    const double Ib_unscaled = opts.getBunchCurrent();
    const double Qb = Ib_unscaled/f_rev;
    const double Ib_scaled = Ib_unscaled/isoscale;
    const double Fk = opts.getStartDistParam();
    const double Iz = opts.getStartDistZoom();

    const unsigned int steps = std::max(opts.getSteps(),1u);
    const unsigned int outstep = opts.getOutSteps();
    const float rotations = opts.getNRotations();
    const double t_d = isoscale*opts.getDampingTime();
    const double dt = 1.0/(fs*steps);
    const double revolutionpart = f0*dt;
    const double t_sync_unscaled = 1.0/fs_unscaled;

    /* angle of one rotation step (in rad)
     * (angle = 2*pi corresponds to 1 synchrotron period)
     */
    const meshaxis_t angle = 2*M_PI/steps;

    std::string startdistfile = opts.getStartDistFile();
    double shield = 0;
    double Ith = 0;
    double S_csr = 0;

    if (gap!=0) {
        if (gap>0) {
            shield = bl*std::sqrt(R_bend)*std::pow(gap,-3./2.);
        }

        const double Inorm = physcons::IAlfven/physcons::me*2*M_PI
                           * std::pow(dE*fs/f0,2)/V/H* std::pow(bl/R_bend,1./3.);

        Ith = Inorm * (0.5+0.34*shield);

        S_csr = Ib_scaled/Inorm;

        if (verbose) {
            sstream.str("");
            sstream << std::fixed << shield;
            Display::printText("Shielding parameter (g=gap):   "
                               +sstream.str());
            if (gap>0) {
                shield = bl*std::sqrt(R_bend)*std::pow(gap/2,-3./2.);
            }
            sstream.str("");
            sstream << std::fixed << shield;
            Display::printText("Shielding parameter (h=height/2): "
                               +sstream.str());
            sstream.str("");
            sstream << std::fixed << S_csr;
            if (Ib_scaled > Ith) {
                sstream << " (> " << 0.5+0.12*shield << ')';
            } else {
                sstream << " (< " << 0.5+0.12*shield << ')';
            }
            Display::printText("CSR strength: "
                               +sstream.str());
            sstream.str("");
            sstream << std::scientific << Ith*isoscale;
            Display::printText("BBT (scaling-law) threshold current at "
                               +sstream.str()+" A.");
            sstream.str("");
            sstream << std::scientific << 1/t_d/fs/(2*M_PI);
            Display::printText("Damping beta: " +sstream.str());
        }
    }

    if (verbose) {
        sstream.str("");
        sstream << std::fixed << 1/revolutionpart;
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
        sstream.precision(2);
        sstream << std::fixed << Fk << ", zoom=" << Iz;
        Display::printText("Generating initial distribution with F(k)="
                           +sstream.str()+".");
        mesh1 = new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,
                               Qb,Ib_unscaled,bl,dE,Fk,Iz);
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
        } catch ( const png::error &e ) {
            std::cerr << "Problem loading " << startdistfile
                      << ": " << e.what() << std::endl;
            return EXIT_SUCCESS;
        } catch (...) {
            std::cerr << "Error loading initial distribution from \""
                      << startdistfile << "\".";
            return EXIT_FAILURE;
        }

        if (image.get_width() == image.get_height()) {
            if (ps_size != image.get_width()) {
                std::cerr << startdistfile
                          << " does not match set GridSize." << std::endl;

                return EXIT_SUCCESS;
            }

            mesh1 = new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,
                                   Qb,Ib_unscaled,bl);

            for (unsigned int x=0; x<ps_size; x++) {
                for (unsigned int y=0; y<ps_size; y++) {
                    (*mesh1)[x][y] = image[ps_size-y-1][x]/float(UINT16_MAX);
                }
            }
            // normalize integral to 1
            mesh1->normalize();

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
    if (  isOfFileType(".h5",startdistfile)
       || isOfFileType(".hdf5",startdistfile) ) {
        try {
            mesh1 = new PhaseSpace(HDF5File::readPhaseSpace(startdistfile,
                                                            qmin,qmax,
                                                            pmin,pmax,
                                                            Qb,Ib_unscaled,
                                                            bl,dE));
        } catch (const std::exception& ex) {
            std::cerr << "Error loading initial distribution from \""
                      << startdistfile << "\":"
                      << ex.what() << std::endl;
            return EXIT_SUCCESS;
        } catch (const H5::Exception& ex) {
            ex.printErrorStack();
            return EXIT_SUCCESS;
        } catch (...) {
            std::cerr << "Error loading initial distribution from \""
                      << startdistfile << "\".";
            return EXIT_FAILURE;
        }

        if (ps_size != mesh1->nMeshCells(0)) {
            std::cerr << startdistfile
                      << " does not match set GridSize." << std::endl;

            delete mesh1;
            return EXIT_SUCCESS;
        }
        mesh1->syncCLMem(clCopyDirection::cpu2dev);
    } else
    #endif
    if (isOfFileType(".txt",startdistfile)) {
        ps_size = opts.getMeshSize();
        mesh1 = new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,Qb,Ib_unscaled,bl);

        std::ifstream ifs;
        try {
            ifs.open(startdistfile);
        } catch (const std::exception& ex) {
            std::cerr << "Error loading initial distribution from \""
                      << startdistfile << "\":"
                      << ex.what() << std::endl;
            return EXIT_SUCCESS;
        } catch (...) {
            std::cerr << "Error loading initial distribution from \""
                      << startdistfile << "\".";
            return EXIT_FAILURE;
        }

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

        // normalize integral to 1
        mesh1->normalize();
        mesh1->syncCLMem(clCopyDirection::cpu2dev);
    } else {
        Display::printText("Unknown format of input file. Will now quit.");
        return EXIT_SUCCESS;
    }
    }

    // find highest peak (for display)
    meshdata_t maxval = std::numeric_limits<meshdata_t>::min();
    for (unsigned int x=0; x<ps_size; x++) {
        for (unsigned int y=0; y<ps_size; y++) {
            maxval = std::max(maxval,(*mesh1)[x][y]);
        }
    }


    #ifdef INOVESA_USE_GUI
    Plot2DLine* bpv = nullptr;
    Plot3DColormap* psv = nullptr;
    Plot2DLine* wpv = nullptr;

    std::vector<float> csrlog;
    Plot2DLine* history = nullptr;
    if (gui) {
        try {
            psv = new Plot3DColormap(maxval);
            display->addElement(psv);
            psv->createTexture(mesh1);
            display->draw();
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(psv);
            delete psv;
            psv = nullptr;
            gui = false;
        }
    }
    #endif // INOVESSA_USE_GUI

    if (verbose) {
    sstream.str("");
    sstream << std::scientific << maxval*Ib_scaled/f0/physcons::e;
    Display::printText("Maximum particles per grid cell is "
                       +sstream.str()+".");
    }

    const double padding =std::max(opts.getPadding(),1.0);

    Impedance* impedance = nullptr;
    if (opts.getImpedanceFile() == "") {
        if (gap>0) {
            Display::printText("Will use parallel plates CSR impedance.");
            impedance = new ParallelPlatesCSR(ps_size*padding,f0,
                                     ps_size*vfps::physcons::c/(2*qmax*bl),gap);
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
        rm1 = new RotationMap(mesh2,mesh3,ps_size,ps_size,angle,
                             interpolationtype,interpol_clamp,
                             RotationMap::RotationCoordinates::norm_pm1,
                             ps_size*ps_size);
        break;
    case 2:
    default:
        Display::printText("Building RFKickMap.");
        rm1 = new RFKickMap(mesh2,mesh1,ps_size,ps_size,angle,
                            interpolationtype,interpol_clamp);

        Display::printText("Building DriftMap.");
        rm2 = new DriftMap(mesh1,mesh3,ps_size,ps_size,
                           {{angle,alpha1/alpha0*angle,alpha2/alpha0*angle}},
                           E0,interpolationtype,interpol_clamp);
        break;
    }

    double e1;
    if (t_d > 0) {
        e1 = 2.0/(fs*t_d*steps);
    } else {
        e1=0;
    }

    SourceMap* fpm;
    if (e1 > 0) {
        Display::printText("Building FokkerPlanckMap.");
        fpm = new FokkerPlanckMap( mesh3,mesh1,ps_size,ps_size,
                                   FokkerPlanckMap::FPType::full,e1,
                                   derivationtype);
    } else {
        fpm = new Identity(mesh3,mesh1,ps_size,ps_size);
    }

    ElectricField* field = nullptr;
    SourceMap* wm = nullptr;
    WakeKickMap* wkm = nullptr;
    WakeFunctionMap* wfm = nullptr;
    std::string wakefile = opts.getWakeFile();
    if (wakefile.size() > 4) {
        field = new ElectricField(mesh1,impedance,revolutionpart);
        Display::printText("Reading WakeFunction from "+wakefile+".");
        wfm = new WakeFunctionMap(mesh1,mesh2,ps_size,ps_size,
                                  wakefile,E0,sE,Ib_scaled,dt,
                                  interpolationtype,interpol_clamp);
        wkm = wfm;
    } else {
        Display::printText("Calculating WakePotential.");
        field = new ElectricField(mesh1,impedance,revolutionpart,
                                  Ib_scaled,E0,sE,dt);
        if (gap != 0) {
            Display::printText("Building WakeKickMap.");
            wkm = new WakePotentialMap(mesh1,mesh2,ps_size,ps_size,field,
                                       interpolationtype,interpol_clamp);
        }
    }
    if (wkm != nullptr) {
        wm = wkm;
    } else {
        wm = new Identity(mesh1,mesh2,ps_size,ps_size);
    }

    // load coordinates for particle tracking
    std::vector<PhaseSpace::Position> trackme;
    if (  opts.getParticleTracking() != ""
       && opts.getParticleTracking() != "/dev/null" ) {
        try {
            std::ifstream trackingfile(opts.getParticleTracking());
            meshaxis_t x,y;
            while (trackingfile >> x >> y) {
                trackme.push_back({x,y});
            }
        } catch (std::exception& e) {
            std::cerr << e.what();
            Display::printText("Will not do particle tracking.");
            trackme.clear();
        }
        std::stringstream npart;
        npart << trackme.size();
        Display::printText( "Will do particle tracking with "
                          + npart.str()
                          + " particles.");
    }

    #ifdef INOVESA_USE_GUI
    if (gui) {
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
        try {
            history = new Plot2DLine(std::array<float,3>{{0,0,0}});
            display->addElement(history);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(wpv);
            delete wpv;
            wpv = nullptr;
        }
    }
    #endif // INOVESA_USE_GUI

    #ifdef INOVESA_USE_HDF5
    HDF5File* hdf_file = nullptr;
    if ( isOfFileType(".h5",ofname)
      || isOfFileType(".hdf5",ofname) ) {
        opts.save(ofname+".cfg");
        Display::printText("Saved configuiration to \""+ofname+".cfg\".");
        hdf_file = new HDF5File(ofname,mesh1,field,impedance,wfm,trackme.size(),
                                t_sync_unscaled);
        Display::printText("Will save results to \""+ofname+"\".");
        opts.save(hdf_file);
        hdf_file->addParameterToGroup("/Info","CSRStrength",
                                      H5::PredType::IEEE_F64LE,&S_csr);
        hdf_file->addParameterToGroup("/Info","ShieldingParameter",
                                      H5::PredType::IEEE_F64LE,&shield);
    } else
    #endif // INOVESA_USE_HDF5
    #ifdef INOVESA_USE_PNG
    if ( isOfFileType(".png",ofname)) {
        opts.save(ofname+".cfg");
        Display::printText("Saved configuiration to \""+ofname+".cfg\".");
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

    #ifdef INOVESA_USE_HDF5
    const HDF5File::AppendType h5save =
        opts.getSavePhaseSpace()? HDF5File::AppendType::All:
                                  HDF5File::AppendType::Defaults;
    if (hdf_file != nullptr && h5save == HDF5File::AppendType::Defaults) {
        // save initial phase space
        hdf_file->append(mesh1,HDF5File::AppendType::PhaseSpace);
    }
    #endif

    if (gui) {
        csrlog.resize(std::floor(steps*rotations/outstep)+1,0);
    }

    Display::printText("Starting the simulation.");
    for (unsigned int i=0, outstepnr=0;i<steps*rotations;i++) {
        if (wkm != nullptr) {
            wkm->update();
        }
        if (outstep > 0 && i%outstep == 0) {
            outstepnr++;
            integral_t meshintegral; // normalized charge (should be 1)
            if (renormalize) {
                // works on XProjection (and recalculates it)
                meshintegral = mesh1->normalize();
            } else {
                // works on XProjection (and recalculates it)
                meshintegral = mesh1->integral();
            }
            mesh1->variance(0);
            mesh1->updateYProjection();
            mesh1->variance(1);
            #ifdef INOVESA_USE_CL
            if (OCLH::active) {
                mesh1->syncCLMem(clCopyDirection::dev2cpu);
                if (wkm != nullptr) {
                    wkm->syncCLMem(clCopyDirection::dev2cpu);
                }
            }
            #endif // INOVESA_USE_CL
            #ifdef INOVESA_USE_HDF5
            if (hdf_file != nullptr) {
                hdf_file->appendTime(static_cast<double>(i)
                                /static_cast<double>(steps));
                hdf_file->append(mesh1,h5save);
                field->updateCSR(fc);
                hdf_file->append(field);
                if (wkm != nullptr) {
                    hdf_file->append(wkm);
                }
                hdf_file->append(trackme.data());
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
                if (history != nullptr) {
                    #ifdef INOVESA_USE_HDF5
                    if (hdf_file == nullptr)
                    #endif // INOVESA_USE_HDF5
                    {
                        field->updateCSR(fc);
                    }
                    csrlog[outstepnr] = field->getCSRPower();
                    history->updateLine(csrlog.size(),csrlog.data(),true);
                }
                display->draw();
                if (psv != nullptr) {
                    psv->delTexture();
                }
            }
            #endif // INOVESSA_USE_GUI
            std::stringstream status;
            status.precision(5);
            status << std::setw(6) << static_cast<float>(i)/steps
                   << '/' << rotations;
            status << "\t1-Q/Q_0=" << 1.0 - meshintegral;
            Display::printText(status.str(),2.0f);
        }
        wm->apply();
        wm->applyTo(trackme);
        rm1->apply();
        rm1->applyTo(trackme);
        if (rm2 != nullptr) {
            rm2->apply();
            rm2->applyTo(trackme);
        }
        fpm->apply();
        fpm->applyTo(trackme);
    }

    #ifdef INOVESA_USE_HDF5
    // save final result
    if (hdf_file != nullptr) {
        if (wkm != nullptr) {
            wkm->update();
        }
        mesh1->integral(); // works on XProjection (and recalculates it)
        mesh1->variance(0);
        mesh1->updateYProjection();
        mesh1->variance(1);
        #ifdef INOVESA_USE_CL
        if (OCLH::active) {
            mesh1->syncCLMem(clCopyDirection::dev2cpu);
            if (wkm != nullptr) {
                wkm->syncCLMem(clCopyDirection::dev2cpu);
            }
        }
        #endif // INOVESA_USE_CL
        hdf_file->appendTime(rotations);
        hdf_file->append(mesh1,HDF5File::AppendType::All);
        field->updateCSR(fc);
        hdf_file->append(field);
        if (wkm != nullptr) {
            hdf_file->append(wkm);
        }
        hdf_file->append(trackme.data());
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
    status.precision(5);
    status << std::setw(6) << rotations << '/' << rotations;
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

    Display::printText("Finished.");

    return EXIT_SUCCESS;
}

