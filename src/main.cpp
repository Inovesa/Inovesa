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
#include "IO/GUI/Plot1DLine.hpp"
#include "IO/GUI/Plot2DPoints.hpp"
#include "IO/GUI/Plot3DColormap.hpp"
#include "PS/PhaseSpace.hpp"
#include "PS/PhaseSpaceFactory.hpp"
#include "Z/ImpedanceFactory.hpp"
#include "CL/OpenCLHandler.hpp"
#include "SM/FokkerPlanckMap.hpp"
#include "SM/Identity.hpp"
#include "SM/KickMap.hpp"
#include "SM/DriftMap.hpp"
#include "SM/RFKickMap.hpp"
#include "SM/DynamicRFKickMap.hpp"
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

#include <boost/math/special_functions/sign.hpp>
#include <boost/math/constants/constants.hpp>
using boost::math::constants::two_pi;

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
     * Will count steps in the main simulation loop,
     * but can be used by time dependent variables.
     */
    uint32_t simulationstep = 0;

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


    oclhptr_t oclh;
    #ifdef INOVESA_USE_OPENCL
    if (cldev < 0) {
        OCLH::listCLDevices();
        return EXIT_SUCCESS;
    }
    if (cldev > 0) {
        try {
            oclh = std::make_shared<OCLH>( opts.getCLDevice()-1
                                         #ifdef INOVESA_USE_OPENGL
                                         , opts.showPhaseSpace()
                                         #endif // INOVESA_USE_OPENGL
                                         );
        } catch (cl::Error& e) {
            Display::printText(e.what());
            Display::printText("Will fall back to sequential version.");
            oclh.reset();
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

    // energy of reference particle (in eV)
    const auto E0 = opts.getBeamEnergy();

    // absolute energy spread (in eV)
    const auto dE = sE*E0;

    /*
     * Only positive benging radii will be used.
     * Otherwise radius will be deduced from revolution frequency
     * (based on iso-magnetic ring model).
     */
    const auto use_set_bend = (opts.getBendingRadius()>0);

    // revolution frequency (in Hz)
    const auto f_rev = opts.getRevolutionFrequency();

    const auto R_bend = use_set_bend
            ? opts.getBendingRadius()
            : physcons::c/(two_pi<double>()*f_rev);

    const frequency_t fc = opts.getCutoffFrequency();
    const auto harmonic_number = opts.getHarmonicNumber();
    const auto f_RF = f_rev*harmonic_number;
    const double gap = opts.getVacuumChamberGap();
    const double V_RF = opts.getRFVoltage();

    const auto lorentzgamma = E0/physcons::me;
    const auto V0 = physcons::e*std::pow(lorentzgamma,4)
                  / (3* physcons::epsilon0*R_bend);

    const auto W0 = V0*physcons::e;

    const double V_eff = std::sqrt(V_RF*V_RF-V0*V0);

    double fs = opts.getSyncFreq();
    meshaxis_t alpha0_tmp = opts.getAlpha0();

    // non-zero f_s will be used, zero implies usage of alpha0
    if (fs == 0) {
        fs = f_rev*std::sqrt(alpha0_tmp*harmonic_number*V_RF/(two_pi<double>()*E0));
    } else {
        alpha0_tmp = boost::math::sign(fs)*two_pi<double>()*E0
                   / (harmonic_number*V_RF)*std::pow(fs/f_rev,2);
    }

    const meshaxis_t alpha0 = alpha0_tmp;
    const meshaxis_t alpha1 = opts.getAlpha1();
    const meshaxis_t alpha2 = opts.getAlpha2();


    // natural RMS bunch length
    const double bl = physcons::c*dE/harmonic_number/std::pow(f_rev,2.0)/V_eff*fs;

    const double Ib = opts.getBunchCurrent();
    const double Qb = Ib/f_rev;
    const double Iz = opts.getStartDistZoom();

    const double steps = (opts.getStepsPerTrev()>0)
            ? opts.getStepsPerTrev()*f_rev/fs
            : std::max(opts.getStepsPerTsync(),1u);
    const auto outstep = opts.getOutSteps();
    const float rotations = opts.getNRotations();

    const auto calc_damp = E0*physcons::e/W0/f_rev;

    const auto set_damp = opts.getDampingTime();

    const auto t_damp = (set_damp < 0)? calc_damp : set_damp;

    const double dt = 1.0/(fs*steps);
    const double revolutionpart = f_rev*dt;
    const double t_sync = 1.0/fs;

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

    // RF Phase Noise Amplitude
    const auto rf_noise_add = std::max(0.0,
                                opts.getRFPhaseSpread()
                                / 360.0*two_pi<double>()
                                * std::sqrt(revolutionpart)*V_RF
                                / dE*ps_size/pqsize);

    // RF Amplitude Noise
    const auto rf_noise_mul = std::max(0.0,
                                opts.getRFAmplitudeSpread()
                                * std::sqrt(revolutionpart)/2.0);


    // RF phase modulation amplitude
    const auto rf_mod_ampl = std::max(0.0,
                                      opts.getRFPhaseModAmplitude()
                                      /360.0*two_pi<double>()
                                      * std::sqrt(revolutionpart));

    // "time step" for RF phase modulation
    const auto rf_mod_step = opts.getRFPhaseModFrequency()*dt;


    /*
     * angle of one rotation step (in rad)
     * (angle = 2*pi corresponds to 1 synchrotron period)
     */
    const meshaxis_t angle = two_pi<double>()/steps;

    uint32_t laststep=std::ceil(steps*rotations);

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

        const double Inorm = physcons::IAlfven/physcons::me*two_pi<double>()
                           * std::pow(dE*fs/f_rev,2)/V_eff/harmonic_number
                           * std::pow(bl/R_bend,1./3.);

        Ith = Inorm * (0.5+0.34*shield);

        S_csr = Ib/Inorm;

        if (verbose && use_csr) {
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
            if (Ib > Ith) {
                sstream << " (> " << 0.5+0.12*shield << ')';
            } else {
                sstream << " (< " << 0.5+0.12*shield << ')';
            }
            Display::printText("CSR strength: "
                               +sstream.str());
            sstream.str("");
            sstream << std::scientific << Ith;
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
        sstream << std::scientific << fs;
        Display::printText("Synchrotron Frequency: " +sstream.str()+ " Hz");

        if (opts.getStepsPerTrev() == 0) {
            sstream.str("");
            sstream << std::fixed << 1/revolutionpart;
            Display::printText("Doing " +sstream.str()+
                               " simulation steps per revolution period.");
        } else {
            sstream.str("");
            sstream << std::fixed << f_rev/fs/revolutionpart;
            Display::printText("Doing " +sstream.str()+
                               " simulation steps per synchrotron period.");
        }

        sstream.str("");
        auto syncphase = std::asin(V0/V_RF)/two_pi<double>()*360;
        sstream << std::fixed << syncphase;
        Display::printText("Synchronous phase is at "+sstream.str()
                           +" degree (should be small).");

        sstream.str("");
        sstream << std::scientific << calc_damp << " s";
        if (set_damp >= 0) {
            sstream << " (set value: "
                    << std::scientific << set_damp  << " s)";
        }
        Display::printText("Damping time calculated from ring parameters is "
                           +sstream.str() + ".");

        sstream.str("");
        sstream << std::scientific << 1/t_damp/fs/(two_pi<double>());
        Display::printText("Damping beta: " +sstream.str());

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
        grid_t1.reset(new PhaseSpace( ps_size,qmin,qmax,pmin,pmax
                                    , oclh
                                    , Qb,Ib,bl,dE,Iz));
    } else {
        Display::printText("Reading in initial distribution from: \""
                           +startdistfile+'\"');
        #ifdef INOVESA_USE_PNG
        // check for file ending .png
        if (isOfFileType(".png",startdistfile)) {
            grid_t1 = makePSFromPNG( startdistfile,qmin,qmax,pmin,pmax
                                   , oclh, Qb,Ib,bl,dE);
        } else
        #endif // INOVESA_USE_PNG
        #ifdef INOVESA_USE_HDF5
        if (  isOfFileType(".h5",startdistfile)
           || isOfFileType(".hdf5",startdistfile) ) {
            grid_t1 = makePSFromHDF5( startdistfile,opts.getStartDistStep()
                                    , qmin,qmax,pmin,pmax
                                    , oclh
                                    , Qb,Ib,bl,dE);

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
            grid_t1 = makePSFromTXT( startdistfile,opts.getGridSize()
                                   , qmin,qmax,pmin,pmax
                                   , oclh
                                   , Qb,Ib,bl,dE);
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
        sstream << std::scientific << maxval*Ib/f_rev/physcons::e;
        Display::printText("Maximum particles per grid cell is "
                           +sstream.str()+".");
    }

    #ifdef INOVESA_USE_OPENGL
    /**************************************************************************
     *  stuff for (graphical) display                                         *
     **************************************************************************/

    // plot of bunch profile
    std::shared_ptr<Plot1DLine> bpv;

    // plot of particle positions
    std::shared_ptr<Plot2DPoints> ppv;

    // plot of phase space
    std::shared_ptr<Plot3DColormap> psv;

    // plot of wake potential
    std::shared_ptr<Plot1DLine> wpv;

    // plot of CSR power (over time)
    std::vector<float> csrlog(std::ceil(steps*rotations/outstep)+1,0);
    std::shared_ptr<Plot1DLine> history;

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
    if ( rf_noise_add != 0 || rf_noise_mul != 0
      || (rf_mod_ampl != 0 && rf_mod_step != 0)) {
        Display::printText("Building dynamic RFKickMap...");
        sstream.str("");
        sstream << opts.getRFPhaseSpread()/360.0/f_rev/harmonic_number;
        if (rf_noise_add != 0) {
            Display::printText("...including phase noise (spread: "
                               + sstream.str()+" s)");
        }
        if (rf_mod_ampl != 0 && rf_mod_step != 0) {
            sstream.str("");
            sstream << "...including phase modulation ("
                    << opts.getRFPhaseModFrequency()/fs
                    << " fs) of +/-"
                    << opts.getRFPhaseModAmplitude()/360.0/f_rev/harmonic_number
                    << " s";
            Display::printText(sstream.str());
        }
        rm1.reset(new DynamicRFKickMap( grid_t2, grid_t1,ps_size, ps_size
                                      , revolutionpart, V_eff, f_RF, V0
                                      , rf_noise_add, rf_noise_mul
                                      , rf_mod_ampl,rf_mod_step
                                      , &simulationstep, laststep
                                      , interpolationtype,interpol_clamp
                                      , oclh
                                      ));
    } else {
        Display::printText("Building static RFKickMap.");
        rm1.reset(new RFKickMap( grid_t2,grid_t1,ps_size,ps_size
                               , revolutionpart, V_eff, f_RF, V0
                               , interpolationtype,interpol_clamp
                               , oclh
                               ));
    }

    Display::printText("Building DriftMap.");
    rm2.reset(new DriftMap( grid_t1,grid_t3,ps_size,ps_size
                          , {{angle,alpha1/alpha0*angle,alpha2/alpha0*angle}}
                          , E0,interpolationtype,interpol_clamp
                          , oclh
                          ));

    // time constant for damping and diffusion
    const timeaxis_t  e1 = (t_damp > 0) ? 2.0/(fs*t_damp*steps) : 0;

    // SourceMap for damping and diffusion
    SourceMap* fpm;
    if (e1 > 0) {
        Display::printText("Building FokkerPlanckMap.");
        fpm = new FokkerPlanckMap( grid_t3,grid_t1,ps_size,ps_size
                                 , FokkerPlanckMap::FPType::full,e1
                                 , derivationtype
                                 , oclh
                                 );
    } else {
        fpm = new Identity(grid_t3,grid_t1,ps_size,ps_size, oclh);
    }



    /*
     * Note: There are two used impedances,
     * one for beam dynamics and one for CSR.
     */

    Display::printText("For beam dynamics computation:");
    std::shared_ptr<Impedance> wake_impedance
            = vfps::makeImpedance( nfreqs
                                 , oclh
                                 , fmax,R_bend,f_rev,gap,use_csr
                                 , s,xi,collimator_radius,impedance_file);

    Display::printText("For CSR computation:");
    std::shared_ptr<Impedance> rdtn_impedance
            = vfps::makeImpedance( nfreqs
                                 , oclh
                                 , fmax,R_bend,f_rev,(gap>0)?gap:-1);


    // field for radiation (not for self-interaction)
    ElectricField rdtn_field( grid_t1,rdtn_impedance
                            , oclh
                            , f_rev,revolutionpart);

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
        wfm = new WakeFunctionMap( grid_t1,grid_t2,ps_size,ps_size
                                 , wakefile,E0,sE,Ib,dt
                                 , interpolationtype,interpol_clamp
                                 , oclh
                                 );
        wkm = wfm;
    } else {
        if (wake_impedance != nullptr) {
            Display::printText("Calculating WakePotential.");
            wake_field = new ElectricField( grid_t1,wake_impedance
                                          , oclh
                                          ,f_rev
                                          , revolutionpart, Ib,E0,sE,dt
                                          );

            Display::printText("Building WakeKickMap.");
            wkm = new WakePotentialMap( grid_t1,grid_t2,ps_size,ps_size
                                      , wake_field ,interpolationtype
                                      , interpol_clamp
                                      , oclh
                                      );
        }
    }
    if (wkm != nullptr) {
        wm = wkm;
    } else {
        wm = new Identity( grid_t1,grid_t2,ps_size,ps_size,oclh);
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
            bpv.reset(new Plot1DLine( std::array<float,3>{{1,0,0}},ps_size
                                    , Plot1DLine::orientation::horizontal));
            display->addElement(bpv);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            display->takeElement(bpv);
            bpv.reset();
        }
        if (!trackme.empty()) {
            try {
                ppv.reset(new Plot2DPoints(std::array<float,3>{{1,1,1}},
                                           ps_size,ps_size));
                display->addElement(ppv);
            } catch (std::exception &e) {
                std::cerr << e.what() << std::endl;
                display->takeElement(ppv);
                ppv.reset();
            }
        }
        if (wkm != nullptr) {
            try {
                wpv.reset(new Plot1DLine( std::array<float,3>{{0,0,1}},ps_size
                                        , Plot1DLine::orientation::horizontal));
                display->addElement(wpv);
            } catch (std::exception &e) {
                std::cerr << e.what() << std::endl;
                display->takeElement(wpv);
                wpv.reset();
            }
        }
        try {
            history.reset(new Plot1DLine( std::array<float,3>{{0,0,0}}
                                        , csrlog.size()
                                        , Plot1DLine::orientation::vertical));
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
                                wfm,trackme.size(), t_sync,f_rev);
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
    const auto h5save = opts.getSavePhaseSpace();
    // end of preparation to save results


    if (hdf_file != nullptr && h5save == 0) {
        // save initial phase space (if not saved anyways)
        hdf_file->append(*grid_t1,0,HDF5File::AppendType::PhaseSpace);
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
    if (oclh) {
        oclh->finish();
    }
    #endif // INOVESA_USE_OPENCL

    /*
     * main simulation loop
     * (everything inside this loop will be run a multitude of times)
     */
    uint32_t outstepnr=0;
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

            // works on XProjection
            grid_t1->getIntegral();
            grid_t1->variance(0);
            grid_t1->updateYProjection();
            grid_t1->variance(1);
            #ifdef INOVESA_USE_OPENCL
            if (oclh) {
                grid_t1->syncCLMem(clCopyDirection::dev2cpu);
                if (wkm != nullptr) {
                    wkm->syncCLMem(clCopyDirection::dev2cpu);
                }
            }
            #endif // INOVESA_USE_OPENCL
            #ifdef INOVESA_USE_HDF5
            if (hdf_file != nullptr) {
                HDF5File::AppendType at =
                        (h5save > 0 && outstepnr%h5save == 0)
                        ? HDF5File::AppendType::All
                        : HDF5File::AppendType::Defaults;


                hdf_file->append(*grid_t1,
                        static_cast<double>(simulationstep)/steps, at);
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
                    rm2->applyTo(allpos);
                    fpm->applyTo(allpos);
                    hdf_file->appendSourceMap(allpos.data());
                }
            }
            outstepnr++;
            #endif // INOVESA_USE_HDF5
            #ifdef INOVESA_USE_OPENGL
            if (display != nullptr) {
                if (psv != nullptr) {
                    psv->createTexture(grid_t1);
                }
                if (bpv != nullptr) {
                    bpv->update(grid_t1->getProjection(0));
                }
                if (ppv != nullptr) {
                    ppv->update(trackme);
                }
                if (wpv != nullptr) {
                    wpv->update(wkm->getForce());
                }
                if (history != nullptr) {
                    #ifdef INOVESA_USE_HDF5
                    if (hdf_file == nullptr)
                    #endif // INOVESA_USE_HDF5
                    {
                        rdtn_field.updateCSR(fc);
                    }
                    csrlog[outstepnr] = rdtn_field.getCSRPower();
                    history->update(csrlog.data());
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
        rm2->apply();
        rm2->applyTo(trackme);
        fpm->apply();
        fpm->applyTo(trackme);

        // udate for next time step
        grid_t1->updateXProjection();

        #ifdef INOVESA_USE_OPENCL
        if (oclh) {
            oclh->flush();
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
        if (renormalize > 0 && simulationstep%renormalize == 0) {
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
        if (oclh) {
            grid_t1->syncCLMem(clCopyDirection::dev2cpu);
            if (wkm != nullptr) {
                wkm->syncCLMem(clCopyDirection::dev2cpu);
            }
        }
        #endif // INOVESA_USE_OPENCL
        // for theresult, everything will be saved
        hdf_file->append(*grid_t1,
                         static_cast<double>(simulationstep)/steps,
                         HDF5File::AppendType::All);
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

