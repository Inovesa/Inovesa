/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
 * Copyright (c) 2018: Patrick Schreiber                                      *
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

#include "defines.hpp"
#include "MessageStrings.hpp"
#include "IO/Display.hpp"
#include "IO/FSPath.hpp"
#include "IO/GUI/Plot2DLine.hpp"
#include "IO/GUI/Plot3DColormap.hpp"
#include "PS/PhaseSpace.hpp"
#include "PS/PhaseSpaceFactory.hpp"
#include "Z/ImpedanceFactory.hpp"
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
#include <memory>
#include <sstream>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

using namespace vfps;

#ifdef INOVESA_ENABLE_INTERRUPT
#include<csignal> // for SIGINT handling

void SIGINT_handler(int) {
    Display::abort = true;
}
#endif // INOVESA_ENABLE_INTERRUPT

int main(int argc, char** argv)
{
    /*
     * Starting time is initialize at the very first moment
     * to have correct timing information, e.g. in the log files.
     * This is a design decision: The program would run as well
     * with a sligtly shifted starting time value.
     */
    Display::start_time = std::chrono::system_clock::now();

    /*
     * Program options might be such that the program does not have
     * to be run at all. As config files (read in based on the command line
     * options) might be errorous, propper error handling is important here.
     */
    ProgramOptions opts;
    try {
        if (!opts.parse(argc,argv)) {
            return EXIT_SUCCESS;
        }
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // see documentation of make_display(...)
    auto cldev = opts.getCLDevice();
    std::string ofname = opts.getOutFile();

    if (ofname.empty() && !opts.getForceRun() && cldev >= 0
        #ifdef INOVESA_USE_OPENGL
        && !opts.showPhaseSpace()
        #endif // INOVESA_USE_OPENGL
       ) {
        std::cout << "Nothing to do. Set at least one of "
                     #ifdef INOVESA_USE_OPENGL
                     " 'gui',"
                     #endif // INOVESA_USE_OPENGL
                     " 'output', or"
                     " 'run_anyway'." << std::endl;
        return EXIT_SUCCESS;
    }

    auto display = (cldev < 0)
                 ? nullptr
                 : make_display( ofname
                               #ifdef INOVESA_USE_OPENGL
                               , opts.showPhaseSpace()
                               , opts.getOpenGLVersion()
                               #endif // INOVESA_USE_OPENGL
                               );

    #ifdef INOVESA_USE_OPENCL
    if (cldev < 0) {
        OCLH::listCLDevices();
        return EXIT_SUCCESS;
    }

    OCLH::active = (cldev > 0);
    if (OCLH::active) {
        try {
            OCLH::prepareCLEnvironment( opts.getCLDevice()-1
                                       #ifdef INOVESA_USE_OPENGL
                                       , opts.showPhaseSpace()
                                       #endif // INOVESA_USE_OPENGL
                                       );
            std::atexit(OCLH::teardownCLEnvironment);
        } catch (cl::Error& e) {
            Display::printText(e.what());
            Display::printText("Will fall back to sequential version.");
            OCLH::active = false;
        }
    }
    #endif // INOVESA_USE_OPENCL

    // here follow a lot of settings and options

    const auto derivationtype = static_cast<FokkerPlanckMap::DerivationType>
            (opts.getDerivationType());

    const auto interpolationtype = static_cast<SourceMap::InterpolationType>
            (opts.getInterpolationPoints());

    const bool interpol_clamp = opts.getInterpolationClamped();
    const bool verbose = opts.getVerbosity();
    const auto renormalize = opts.getRenormalizeCharge();

    const auto save_sourcemap = opts.getSaveSourceMap();

    const meshindex_t ps_size = opts.getGridSize();
    const double pqsize = opts.getPhaseSpaceSize();
    const double qcenter = -opts.getPSShiftX()*pqsize/(ps_size-1);
    const double pcenter = -opts.getPSShiftY()*pqsize/(ps_size-1);
    const double pqhalf = pqsize/2;
    const double qmax = qcenter + pqhalf;
    const double qmin = qcenter - pqhalf;
    const double pmax = pcenter + pqhalf;
    const double pmin = pcenter - pqhalf;

    // relative energy spread
    const auto sE = opts.getEnergySpread();

    // energy of reference particle
    const auto E0 = opts.getBeamEnergy();

    // absolute energy spread
    const auto dE = sE*E0;

    // revolution frequency
    const auto f_rev = opts.getRevolutionFrequency();

    /*
     * Only positive benging radii will be used.
     * Otherwise radius will be deduced from revolution frequency
     * (based on iso-magnetic ring model).
     */
    const auto use_set_bend = (opts.getBendingRadius()>0);
    const auto R_bend = use_set_bend
            ? opts.getBendingRadius()
            : physcons::c/(2*pi<double>()*f_rev);

    /*
     * Revolution frequency an iso-magnetic ring
     * with same bending radius would have.
     * It is used as fundamental frequency for impedances.
     */
    const auto f0 = !use_set_bend
            ? f_rev
            : physcons::c/(2*pi<double>()*R_bend);

    // scaling for isomagnetic approximation, defined to be <= 1
    const double isoscale = f_rev/f0;

    const frequency_t fc = opts.getCutoffFrequency();
    const double H_unscaled = opts.getHarmonicNumber();
    const double H = isoscale*H_unscaled;
    const double gap = opts.getVacuumChamberGap();
    const double V = opts.getRFVoltage();

    double fs_tmp = opts.getSyncFreq();
    meshaxis_t alpha0_tmp = opts.getAlpha0();

    // non-zero f_s will be used, zero implies usage of alpha0
    if (fs_tmp == 0) {
        fs_tmp = f_rev*std::sqrt(alpha0_tmp*H_unscaled*V/(2*pi<double>()*E0));
    } else {
        // alpha0 should have same sign as fs
        auto sign = (fs_tmp > 0) ? 1 : -1;

        alpha0_tmp = sign*2*pi<double>()*E0/(H_unscaled*V)*std::pow(fs_tmp/f_rev,2);
    }

    // synchrotron frequency (comparable to real storage ring)
    const double fs_unscaled = fs_tmp;

    // synchrotron frequency (isomagnetic ring)
    const double fs = fs_unscaled/isoscale;

    const meshaxis_t alpha0 = alpha0_tmp;
    const meshaxis_t alpha1 = opts.getAlpha1();
    const meshaxis_t alpha2 = opts.getAlpha2();


    // natural RMS bunch length
    const double bl = physcons::c*dE/H/std::pow(f0,2.0)/V*fs;

    const double Ib_unscaled = opts.getBunchCurrent();
    const double Qb = Ib_unscaled/f_rev;
    const double Ib_scaled = Ib_unscaled/isoscale;
    const double Iz = opts.getStartDistZoom();

    const unsigned int steps = std::max(opts.getSteps(),1u);
    const unsigned int outstep = opts.getOutSteps();
    const float rotations = opts.getNRotations();
    const double t_d = isoscale*opts.getDampingTime();
    const double dt = 1.0/(fs*steps);
    const double revolutionpart = f0*dt;
    const double t_sync_unscaled = 1.0/fs_unscaled;

    const double padding =std::max(opts.getPadding(),1.0);
    const frequency_t fmax = ps_size*vfps::physcons::c/(pqsize*bl);
    const size_t nfreqs = opts.getRoundPadding() ?
                    Impedance::upper_power_of_two(ps_size*padding) :
                    ps_size*padding;
    const auto s = opts.getWallConductivity();
    const auto xi = opts.getWallSusceptibility();
    const auto collimator_radius = opts.getCollimatorRadius();
    const auto impedance_file = opts.getImpedanceFile();
    const auto use_csr = opts.getUseCSR();

    /*
     * angle of one rotation step (in rad)
     * (angle = 2*pi corresponds to 1 synchrotron period)
     */
    const meshaxis_t angle = 2*pi<double>()/steps;

    std::string startdistfile = opts.getStartDistFile();


    /*
     * needed for output text in the main function
     * (move functionality to Display at some point)
     */
    std::stringstream sstream;

    /**************************************************************************
     * Up next: Printing information (dynamics estimation) for upcoming run.  *
     * Only some part has its own context because information will go to      *
     * the results file.                                                      *
     **************************************************************************/
    double shield = 0;
    double Ith = 0;
    double S_csr = 0;

    { // context of information printing, not needed in the program
    if (gap!=0) {
        if (gap>0) {
            shield = bl*std::sqrt(R_bend)*std::pow(gap,-3./2.);
        }

        const double Inorm = physcons::IAlfven/physcons::me*2*pi<double>()
                           * std::pow(dE*fs/f0,2)/V/H
                           * std::pow(bl/R_bend,1./3.);

        Ith = Inorm * (0.5+0.34*shield);

        S_csr = Ib_scaled/Inorm;

        if (verbose) {
            sstream.str("");
            sstream << std::fixed << shield;
            Display::printText("Shielding parameter (g=gap): "
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
        }
    }

    if (verbose) {
        sstream.str("");
        sstream << std::scientific
                << bl << " m ("
                << bl/physcons::c << " s)"
                ;
        Display::printText("Natural bunch length is "
                           +sstream.str());

        sstream.str("");
        sstream << std::scientific << fs_unscaled;
        Display::printText("Synchrotron Frequency: " +sstream.str()+ " Hz");

        sstream.str("");
        sstream << std::scientific << 1/t_d/fs/(2*pi<double>());
        Display::printText("Damping beta: " +sstream.str());

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
    } // end of context of information printing


    /*
     * There are three phase space grids:
     * grid_t1, grid_t2, and grid_t3
     * This will allow to work with odd and even numbers of SourceMaps
     * f(x,y,t) -> f(x,y,t+dt) with fixed source and destination.
     * When memory usage is crytical, it might be worth to change this
     * to two grids (only).
     *
     * This first grid (grid_t1) will be initialized and
     * copied for the other ones.
     */
    std::shared_ptr<PhaseSpace> grid_t1;


    /*
     * initialization of mesh1
     *
     * TODO: It can be considered ugly to do this in main(),
     * so initialization might be moved to a factory function
     * at some point.
     */
    if (startdistfile.empty()) {
        if (ps_size == 0) {
            Display::printText("Please give file for initial distribution "
                               "or size of target mesh > 0.");
        }
        grid_t1.reset(new PhaseSpace(ps_size,qmin,qmax,pmin,pmax,
                               Qb,Ib_unscaled,bl,dE,Iz));
    } else {
        Display::printText("Reading in initial distribution from: \""
                           +startdistfile+'\"');
        #ifdef INOVESA_USE_PNG
        // check for file ending .png
        if (isOfFileType(".png",startdistfile)) {
            grid_t1 = makePSFromPNG(startdistfile,qmin,qmax,pmin,pmax,
                                Qb,Ib_unscaled,bl,dE);
        } else
        #endif // INOVESA_USE_PNG
        #ifdef INOVESA_USE_HDF5
        if (  isOfFileType(".h5",startdistfile)
           || isOfFileType(".hdf5",startdistfile) ) {
            grid_t1 = makePSFromHDF5(startdistfile,opts.getStartDistStep(),
                                   qmin,qmax,pmin,pmax,
                                   Qb,Ib_unscaled,bl,dE);

            if (grid_t1 == nullptr) {
                return EXIT_SUCCESS;
            }

            if (ps_size != grid_t1->nMeshCells(0)) {
                std::cerr << startdistfile
                          << " does not match set GridSize." << std::endl;

                return EXIT_SUCCESS;
            }
        } else
        #endif
        if (isOfFileType(".txt",startdistfile)) {
            grid_t1 = makePSFromTXT(startdistfile,opts.getGridSize(),
                                  qmin,qmax,pmin,pmax,
                                  Qb,Ib_unscaled,bl,dE);
        } else {
            Display::printText("Unknown format of input file. Will now quit.");
            return EXIT_SUCCESS;
        }
    }

    // an initial renormalization might be applied
    if (renormalize >= 0) {
        grid_t1->updateXProjection();

        grid_t1->normalize(); // works on XProjection
    }

    auto grid_t2 = std::make_shared<PhaseSpace>(*grid_t1);
    auto grid_t3 = std::make_shared<PhaseSpace>(*grid_t1);

    // find highest peak for display (and information in the log)
    meshdata_t maxval = std::numeric_limits<meshdata_t>::min();
    for (unsigned int x=0; x<ps_size; x++) {
        for (unsigned int y=0; y<ps_size; y++) {
            maxval = std::max(maxval,(*grid_t1)[x][y]);
        }
    }

    if (verbose) {
        sstream.str("");
        sstream << std::scientific << maxval*Ib_scaled/f0/physcons::e;
        Display::printText("Maximum particles per grid cell is "
                           +sstream.str()+".");
    }

    #ifdef INOVESA_USE_OPENGL
    /**************************************************************************
     *  stuff for (graphical) display                                         *
     **************************************************************************/

    // plot of bunch profile
    std::shared_ptr<Plot2DLine> bpv;

    // plot of phase space
    std::shared_ptr<Plot3DColormap> psv;

    // plot of wake potential
    std::shared_ptr<Plot2DLine> wpv;

    // plot of CSR power (over time)
    std::vector<float> csrlog(std::ceil(steps*rotations/outstep)+1,0);
    std::shared_ptr<Plot2DLine> history;

    /*
     * Plot of phase space is initialized here already,
     * so that users do not have to wait for the actual simulation
     * to initialize: For a first view, the initial phase space is enough.
     */
    if (display != nullptr) {
        try {
            psv.reset(new Plot3DColormap(maxval));
            display->addElement(psv);
            psv->createTexture(grid_t1);
            display->draw();
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(psv);
            psv.reset();
        }
    }
    #endif // INOVESSA_USE_GUI


    // rotation map(s): two, in case of Manhattan rotation
    std::unique_ptr<SourceMap> rm1;
    std::unique_ptr<SourceMap> rm2;
    const uint_fast8_t rotationtype = opts.getRotationType();
    switch (rotationtype) {
    case 0:
        Display::printText("Initializing RotationMap.");
        rm1.reset(new RotationMap(grid_t2,grid_t3,ps_size,ps_size,angle,
                             interpolationtype,interpol_clamp,0));
        break;
    case 1:
        Display::printText("Building RotationMap.");
        rm1.reset(new RotationMap(grid_t2,grid_t3,ps_size,ps_size,angle,
                                  interpolationtype,interpol_clamp,
                                  ps_size*ps_size));
        break;
    case 2:
    default:
        Display::printText("Building RFKickMap.");
        rm1.reset(new RFKickMap(grid_t2,grid_t1,ps_size,ps_size,angle,
                            interpolationtype,interpol_clamp));

        Display::printText("Building DriftMap.");
        rm2.reset(new DriftMap(grid_t1,grid_t3,ps_size,ps_size,
                           {{angle,alpha1/alpha0*angle,alpha2/alpha0*angle}},
                           E0,interpolationtype,interpol_clamp));
        break;
    }
    if (rotationtype != 2 && (alpha1 != 0.0 || alpha2 != 0.0)) {
        Display::printText("Warning: Nonlinear momentum compaction "
                           "incompatible with classical rotation.");
    }

    // time constant for damping and diffusion
    const timeaxis_t  e1 = (t_d > 0) ? 2.0/(fs*t_d*steps) : 0;

    // SourceMap for damping and diffusion
    SourceMap* fpm;
    if (e1 > 0) {
        Display::printText("Building FokkerPlanckMap.");
        fpm = new FokkerPlanckMap( grid_t3,grid_t1,ps_size,ps_size,
                                   FokkerPlanckMap::FPType::full,e1,
                                   derivationtype);
    } else {
        fpm = new Identity(grid_t3,grid_t1,ps_size,ps_size);
    }



    /*
     * Note: There are two used impedances,
     * one for beam dynamics and one for CSR.
     */

    Display::printText("For beam dynamics computation:");
    std::shared_ptr<Impedance> wake_impedance
            = vfps::makeImpedance(nfreqs,fmax,f0,f_rev,gap,use_csr,
                                  s,xi,collimator_radius,impedance_file);

    Display::printText("For CSR computation:");
    std::shared_ptr<Impedance> rdtn_impedance
            = vfps::makeImpedance(nfreqs,fmax,f0,f_rev,(gap>0)?gap:-1);


    // field for radiation (not for self-interaction)
    ElectricField rdtn_field(grid_t1,rdtn_impedance,f_rev,revolutionpart);

    /**************************************************************************
     * Part modeling the self-interaction of the electron-bunch.              *
     **************************************************************************/

    ElectricField* wake_field = nullptr;

    // (generic) source map, will be executed in the main loop
    SourceMap* wm = nullptr;

    // depending if working time or frequency domain, only one might be used
    WakeKickMap* wkm = nullptr;
    WakeFunctionMap* wfm = nullptr;

    std::string wakefile = opts.getWakeFile();
    if (wakefile.size() > 4) {
        Display::printText("Reading WakeFunction from "+wakefile+".");
        wfm = new WakeFunctionMap(grid_t1,grid_t2,ps_size,ps_size,
                                  wakefile,E0,sE,Ib_scaled,dt,
                                  interpolationtype,interpol_clamp);
        wkm = wfm;
    } else {
        if (wake_impedance != nullptr) {
            Display::printText("Calculating WakePotential.");
            wake_field = new ElectricField(grid_t1,wake_impedance,f_rev,
                                           revolutionpart, Ib_scaled,E0,sE,dt);

            Display::printText("Building WakeKickMap.");
            wkm = new WakePotentialMap(grid_t1,grid_t2,ps_size,ps_size,wake_field,
                                     interpolationtype,interpol_clamp);
        }
    }
    if (wkm != nullptr) {
        wm = wkm;
    } else {
        wm = new Identity(grid_t1,grid_t2,ps_size,ps_size);
    }

    /* Load coordinates for particle tracking.
     * Particle tracking is for visualization puproses only,
     * actual beam dynamics may not be perfectly accurate.
     */
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

    // initialze the rest of the display elements
    #ifdef INOVESA_USE_OPENGL
    if (display != nullptr) {
        try {
            bpv.reset(new Plot2DLine(std::array<float,3>{{1,0,0}}));
            display->addElement(bpv);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(bpv);
            bpv.reset();
        }
        if (wkm != nullptr) {
            try {
                wpv.reset(new Plot2DLine(std::array<float,3>{{0,0,1}}));
                display->addElement(wpv);
            } catch (std::exception &e) {
                std::cerr << e.what() << std::endl;
                display->takeElement(wpv);
                wpv.reset();
            }
        }
        try {
            history.reset(new Plot2DLine(std::array<float,3>{{0,0,0}}));
            display->addElement(history);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(history);
            history.reset();
        }
    }
    #endif // INOVESA_USE_OPENGL

    /*
     * preparation to save results
     */
    #ifdef INOVESA_USE_HDF5
    HDF5File* hdf_file = nullptr;
    if ( isOfFileType(".h5",ofname)
      || isOfFileType(".hdf5",ofname) ) {
        opts.save(ofname+".cfg");
        Display::printText("Saved configuiration to \""+ofname+".cfg\".");
        hdf_file = new HDF5File(ofname,grid_t1, &rdtn_field, wake_impedance,
                                wfm,trackme.size(), t_sync_unscaled,f_rev);
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
    if ( ofname.empty() ) {
        Display::printText("Will not save results.");
    } else {
        Display::printText("Unkown filetype for output.");
        return EXIT_SUCCESS;
    }

    #ifdef INOVESA_USE_HDF5
    const HDF5File::AppendType h5save =
        opts.getSavePhaseSpace()? HDF5File::AppendType::All:
                                  HDF5File::AppendType::Defaults;
    // end of preparation to save results


    if (hdf_file != nullptr && h5save == HDF5File::AppendType::Defaults) {
        // save initial phase space (if not saved anyways)
        hdf_file->append(*grid_t1,HDF5File::AppendType::PhaseSpace);
    }
    #endif



    Display::printText("Starting the simulation.");

    // time between two status updates (in seconds)
    const auto updatetime = 2.0f;

    /* We claim that simulation starts now (see log output above).
     * To have the first step always displayed, we do it outside the loop
     * there are two pieces of information needed for this (see below). */

    // 1) the integral
    grid_t1->updateXProjection();
    grid_t1->integrate();

    // 2) the energy spread (variance in Y direction)
    grid_t1->updateYProjection();
    grid_t1->variance(1);
    Display::printText(status_string(grid_t1,0,rotations),false);

    #ifdef INOVESA_ENABLE_INTERRUPT
    //Install signal handler for SIGINT
    signal(SIGINT, SIGINT_handler);
    #endif // INOVESA_ENABLE_INTERRUPT

    #ifdef INOVESA_USE_OPENCL
    if (OCLH::active) {
        OCLH::finish();
    }
    #endif // INOVESA_USE_OPENCL

    /*
     * main simulation loop
     * (everything inside this loop will be run a multitude of times)
     */
    unsigned int simulationstep=0;
    unsigned int outstepnr=0;
    unsigned int laststep=steps*rotations;
    while (simulationstep<laststep && !Display::abort) {
        if (wkm != nullptr) {
            // works on XProjection
            wkm->update();
        }
        if (renormalize > 0 && simulationstep%renormalize == 0) {
            // works on XProjection
            grid_t1->normalize();
        } else {
            // works on XProjection
            grid_t1->integrate();
        }

        if (outstep > 0 && simulationstep%outstep == 0) {
            outstepnr++;

            // works on XProjection
            grid_t1->getIntegral();
            grid_t1->variance(0);
            grid_t1->updateYProjection();
            grid_t1->variance(1);
            #ifdef INOVESA_USE_OPENCL
            if (OCLH::active) {
                grid_t1->syncCLMem(clCopyDirection::dev2cpu);
                if (wkm != nullptr) {
                    wkm->syncCLMem(clCopyDirection::dev2cpu);
                }
            }
            #endif // INOVESA_USE_OPENCL
            #ifdef INOVESA_USE_HDF5
            if (hdf_file != nullptr) {
                hdf_file->appendTime(static_cast<double>(simulationstep)
                                /static_cast<double>(steps));
                hdf_file->append(*grid_t1,h5save);
                rdtn_field.updateCSR(fc);
                hdf_file->append(&rdtn_field);
                if (wkm != nullptr) {
                    hdf_file->append(wkm);
                }
                hdf_file->appendTracks(trackme.data());

                if (save_sourcemap) {
                    std::vector<PhaseSpace::Position> allpos;
                    for (float x=0; x<ps_size; x++) {
                        for (float y=0; y<ps_size; y++) {
                            allpos.push_back({x,y});
                        }
                    }
                    wm->applyTo(allpos);
                    rm1->applyTo(allpos);
                    if (rm2 != nullptr) {
                        rm2->applyTo(allpos);
                    }
                    fpm->applyTo(allpos);
                    hdf_file->appendSourceMap(allpos.data());
                }
            }
            #endif // INOVESA_USE_HDF5
            #ifdef INOVESA_USE_OPENGL
            if (display != nullptr) {
                if (psv != nullptr) {
                    psv->createTexture(grid_t1);
                }
                if (bpv != nullptr) {
                    bpv->updateLine(grid_t1->nMeshCells(0),
                                    grid_t1->getProjection(0));
                }
                if (wpv != nullptr) {
                    wpv->updateLine(grid_t1->nMeshCells(0),
                                    wkm->getForce());
                }
                if (history != nullptr) {
                    #ifdef INOVESA_USE_HDF5
                    if (hdf_file == nullptr)
                    #endif // INOVESA_USE_HDF5
                    {
                        rdtn_field.updateCSR(fc);
                    }
                    csrlog[outstepnr] = rdtn_field.getCSRPower();
                    history->updateLine(csrlog.size(),csrlog.data(),true);
                }
                display->draw();
                if (psv != nullptr) {
                    psv->delTexture();
                }
            }
            #endif // INOVESSA_USE_GUI
            Display::printText(status_string(grid_t1,static_cast<float>(simulationstep)/steps,
                               rotations),false,updatetime);
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

        // udate for next time step
        grid_t1->updateXProjection();

        #ifdef INOVESA_USE_OPENCL
        if (OCLH::active) {
            OCLH::flush();
        }
        #endif // INOVESA_USE_OPENCL

        simulationstep++;
    } // end of main simulation loop

    #ifdef INOVESA_USE_HDF5
    // save final result
    if (hdf_file != nullptr) {
        if (wkm != nullptr) {
            wkm->update();
        }
        /* Without renormalization at this point
         * the last time step might behave slightly different
         * from the ones before.
         */
        if (renormalize > 0) {
            // works on XProjection
            grid_t1->normalize();
        } else {
            // works on XProjection
            grid_t1->integrate();
        }
        grid_t1->variance(0);
        grid_t1->updateYProjection();
        grid_t1->variance(1);
        #ifdef INOVESA_USE_OPENCL
        if (OCLH::active) {
            grid_t1->syncCLMem(clCopyDirection::dev2cpu);
            if (wkm != nullptr) {
                wkm->syncCLMem(clCopyDirection::dev2cpu);
            }
        }
        #endif // INOVESA_USE_OPENCL
        hdf_file->appendTime(static_cast<double>(simulationstep) /static_cast<double>(steps));

        // for the final result, everything will be saved
        hdf_file->append(*grid_t1,HDF5File::AppendType::All);
        rdtn_field.updateCSR(fc);
        hdf_file->append(&rdtn_field);
        if (wkm != nullptr) {
            hdf_file->append(wkm);
        }
        hdf_file->appendTracks(trackme.data());

        if (save_sourcemap) {
            std::vector<PhaseSpace::Position> allpos;
            for (float x=0; x<ps_size; x++) {
                for (float y=0; y<ps_size; y++) {
                    allpos.push_back({x,y});
                }
            }
            wm->applyTo(allpos);
            rm1->applyTo(allpos);
            if (rm2 != nullptr) {
                rm2->applyTo(allpos);
            }
            fpm->applyTo(allpos);
            hdf_file->appendSourceMap(allpos.data());
        }
    }
    #endif // INOVESA_USE_HDF5
    #ifdef INOVESA_USE_PNG
    if ( isOfFileType(".png",ofname)) {
        meshdata_t maxval = std::numeric_limits<meshdata_t>::min();
        meshdata_t* val = grid_t1->getData();
        for (meshindex_t i=0; i<grid_t1->nMeshCells(); i++) {
            maxval = std::max(val[i],maxval);
        }
        png::image< png::gray_pixel_16 > png_file(ps_size, ps_size);
        for (unsigned int x=0; x<ps_size; x++) {
            for (unsigned int y=0; y<ps_size; y++) {
                png_file[ps_size-y-1][x]=
                        static_cast<png::gray_pixel_16>(
                            std::max((*grid_t1)[x][y],meshdata_t(0))
                            /maxval*float(UINT16_MAX));
            }
        }
        png_file.write(ofname);
    }
    #endif

    // Print the last status.
    Display::printText(status_string(grid_t1, static_cast<float>(simulationstep)/steps, rotations));

    delete wake_field;

    delete wm;
    delete fpm;

    // Print Aborted instead of Finished if it was aborted. Also for log file.
    if(Display::abort) {
        Display::printText("Aborted.");
    } else {
        Display::printText("Finished.");
    }

    return EXIT_SUCCESS;
}

