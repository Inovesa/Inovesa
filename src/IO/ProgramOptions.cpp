// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Patrick Schreiber
 * Copyright (c) Johannes Schestag
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include "CL/OpenCLHandler.hpp"

#include "IO/ProgramOptions.hpp"

vfps::ProgramOptions::ProgramOptions() :
    _configfile("default.cfg"),
    I_b({3e-3f}),
    _physopts("Physical Parameters for Simulation"),
    _proginfoopts("Program Information"),
    _programopts_cli("General Program Parameters"),
    _simulopts("Non-Physical Parameters for Simulation"),
    _visibleopts("Possible Parameters")
{
    _proginfoopts.add_options()
        ("help,h", "print help message")
        ("copyright", "print copyright information")
        ("version", "print version string")
        ("buildinfo", "print information on current build")
    ;
    _physopts.add_options()
        ("alpha0", po::value<decltype(alpha0)>(
             &alpha0)->default_value(static_cast<decltype(alpha0)>
                                     (4e-3f), "4e-3"),
            "Linear Momentum compaction factor (1)")
        ("alpha1", po::value<decltype(alpha1)>(
             &alpha1)->default_value(0),
            "Quadratic Momentum compaction factor (1)")
        ("alpha2", po::value<decltype(alpha2)>(
             &alpha2)->default_value(0),
            "Cubic Momentum compaction factor (1)")
        ("SynchrotronFrequency,f", po::value<decltype(f_s)>(
             &f_s)->default_value(0, "(ignore)"),
            "Synchrotron frequency (Hz), "
            "will overwrite alpha0 when set to a value different from 0")
        ("RevolutionFrequency,F",po::value<decltype(f0)>(
             &f0)->default_value(static_cast<decltype(f0)>(9e6), "9e6"),
            "Revolution frequency (Hz)")
        ("DampingTime,d", po::value<decltype(t_d)>(
             &t_d)->default_value(-1),
            "Damping time (s), "
            "<0: calculate based on other parameters")
        ("HarmonicNumber,H", po::value<decltype(H)>(
             &H)->default_value(50),
            "Harmonic Number (1)")
        ("InitialDistFile,i", po::value<decltype(_startdistfile)>(
             &_startdistfile),
            "might be:\n"
            #if INOVESA_USE_HDF5 == 1
            "\tInovesa result file (.hdf5, .h5)\n"
            #endif // INOVESA_USE_HDF5
            #if INOVESA_USE_PNG == 1
            "\tgrayscale png (.png) file\n"
            #endif // INOVESA_USE_PNG
            "\ttext file (.txt) w/ particle coordinates\n"
            "\t'/dev/null' to explicitly state no read-in")
        ("InitialDistStep",po::value<decltype(_startdiststep)>(
             &_startdiststep)->default_value(-1),
            "Select step of HDF5 file for initial distribution")
        ("InitialDistZoom",po::value<decltype(zoom)>(
             &zoom)->default_value(1),
            "Magnification for generation of initial distribution")
        ("BunchCurrent,I", po::value<std::vector<integral_t>>(
             &I_b)->multitoken(),
            "Bunch currents (A)\n"
            " 0: empty bucket (long gap)\n"
            ">0: bucket containing electrons")
        ("BendingRadius,R", po::value<decltype(r_bend)>(
             &r_bend)->default_value(-1),
            "Bending radius of accelerator (m)\n"
            "negative: calculate from RevolutionFrequency")
        ("BeamEnergy,E", po::value<decltype(E_0)>(
             &E_0)->default_value(1.3e9,"1.3e9"),
            "Beam energy (eV)")
        ("BeamEnergySpread,e", po::value<decltype(s_E)>(
             &s_E)->default_value(4.7e-4,"4.7e-4"),
            "Natural energy spread (relative)")
        ("Impedance,Z", po::value<std::string>(
             &_impedancefile),
            "File containing impedance information.")
        ("VacuumGap,G", po::value<decltype(g)>(
             &g)->default_value(0.03,"0.03"),
            "Full height of vacuum chamber (m)\n"
            "<0: free space CSR\n"
            " 0: no CSR\n"
            ">0: parallel plates CSR\n"
            "(|G| is used as size of the beam pipe for other impedances)")
        ("UseCSR", po::value<bool>(
             &use_csr)->default_value(true),
            "Switch to turn off CSR for VacuumGap != 0")
        ("CollimatorRadius", po::value<decltype(collimator)>(
             &collimator)->default_value(0),
            "Radius of collimator opening (m)\n"
            "<=0: no collimator")
        ("WallConductivity", po::value<decltype(s_c)>(
             &s_c)->default_value(0),
            "Conductivity of the vacuum pipe (S/m)\n"
            "<=0: perfect conductor")
        ("WallSusceptibility", po::value<decltype(xi_wall)>(
             &xi_wall)->default_value(0),
            "Magnetic susceptibility of the vacuum pipe (1)\n"
            "<-1: -1")
        ("CutoffFreq", po::value<decltype(f_c)>(
             &f_c)->default_value(static_cast<decltype(f_c)>(23e9),"23e9"),
            "Beamline cutoff frequency (Hz)")
        ("AcceleratingVoltage,V",
            po::value<decltype(V_RF)>(
             &V_RF)->default_value(static_cast<decltype(V_RF)>(1e6),"1e6"),
            "Accelerating Voltage (V) for one revolution")
        ("LinearRF",
            po::value<decltype(linearRF)>(
             &linearRF)->default_value(true),
            "Use linear model for accelerating voltage (experimental)")
        ("RFAmplitudeSpread",
            po::value<decltype(rf_amplitude_spread)>(
             &rf_amplitude_spread)->default_value(0),
            "Relative accelerating voltage amplitude spread per turn")
        ("RFPhaseSpread",
            po::value<decltype(rf_phase_spread)>(
             &rf_phase_spread)->default_value(0),
            "Absolute accelerating voltage phase spread per turn (degree)")
        ("RFPhaseModAmplitude",
            po::value<decltype(rf_phase_mod_amplitude)>(
             &rf_phase_mod_amplitude)->default_value(0),
            "Accelerating voltage phase modulation amplitude (degree)")
        ("RFPhaseModFrequency",
            po::value<decltype(rf_phase_mod_frequency)>(
             &rf_phase_mod_frequency)->default_value(0),
            "Accelerating voltage phase modulation frequency (Hz)")
    ;
    _programopts_file.add_options()
        ("cldev", po::value<decltype(_cldevice)>(
             &_cldevice)->default_value(0),
            "OpenCL device to use\n('-1' lists available devices)")
        ("ForceOpenGLVersion", po::value<decltype(_glversion)>(
             &_glversion)->default_value(2),
            "Force OpenGL version")
        ("gui,g", po::value<decltype(_showphasespace)>(
             &_showphasespace)->default_value(false),
            "Show phase space view")
        ("output,o",
            po::value<decltype(_outfile)>(
             &_outfile),
            "name of file to safe results.")
        ("outstep,n", po::value<decltype(outsteps)>(
             &outsteps)->default_value(100),
            "Save results every n steps.")
        ("SavePhaseSpace",
            po::value<decltype(_savephasespace)>(
             &_savephasespace)->default_value(0),
            "save every n's outstep's phase space to HDF5 file")
        ("tracking",
            po::value<decltype(_trackingfile)>(
             &_trackingfile)->default_value(""),
            "file containing starting positions (grid points)"
            "of particles to be (pseudo-) tracked")
        ("verbose,v", po::value<bool>(
             &_verbose)->default_value(false),
            "print information more detailed")
        ("run_anyway", po::value<bool>(
             &_forcerun)->default_value(false),
            "set to omit consistency check for parameters")
    ;
    _programopts_cli.add_options()
            ("cldev", po::value<decltype(_cldevice)>(
                 &_cldevice)->default_value(0),
        #if INOVESA_USE_OPENCL == 1
            "OpenCL device to use\n('-1' lists available devices)")
        #else // not INOVESA_USE_OPENCL
            "(ignored in this build)")
        #endif // INOVESA_USE_OPENCL
        ("config,c", po::value<decltype(_configfile)>(
             &_configfile),
            "name of a file containing a configuration.")
        ("ForceOpenGLVersion", po::value<decltype(_glversion)>(
             &_glversion)->default_value(2),
            "Force OpenGL version")
        ("gui,g", po::value<decltype(_showphasespace)>(
             &_showphasespace)->default_value(true)->implicit_value(true),
            "Show phase space view")
        ("output,o",
            po::value<std::string>(
             &_outfile),
            "name of file to safe results.")
        ("outstep,n", po::value<decltype(outsteps)>(
             &outsteps)->default_value(100),
            "Save results every n steps.")
        ("SavePhaseSpace",
            po::value<decltype(_savephasespace)>(
             &_savephasespace)->default_value(0),
            "save every n's outstep's phase space to HDF5 file")
        ("tracking",
            po::value<decltype(_trackingfile)>(
             &_trackingfile)->default_value(""),
            "file containing starting positions (grid points)"
            "of particles to be (pseudo-) tracked")
        ("verbose,v", po::value<decltype(_verbose)>(
             &_verbose)->default_value(false)->implicit_value(true),
            "print information more detailed")
        ("run_anyway", po::value<decltype(_forcerun)>(
             &_forcerun)->default_value(false)->implicit_value(true),
            "set to omit consistency check for parameters")
    ;
    _simulopts.add_options()
        ("StepsPerTs,N", po::value<decltype(steps_per_Ts)>(
             &steps_per_Ts)->default_value(1000),
            "Simulation steps for one synchrotron period")
        ("StepsPerRevolution", po::value<decltype(steps_per_Trev)>(
             &steps_per_Trev)->default_value(0),
            "Simulation steps for one revolution (overwrites StepsPerTs)")
        ("padding,p", po::value<decltype(padding)>(
             &padding)->default_value(static_cast<decltype(padding)>(8.0)),
            "Factor for zero padding of single bunch profile(s)")
        ("RoundPadding", po::value<decltype(roundpadding)>(
             &roundpadding)->default_value(true),
            "Always do zero padding up to 2 to the power of N")
        ("PhaseSpaceSize,P", po::value<decltype(pq_size)>(
             &pq_size)->default_value(12),
            "Size of phase space")
        ("PhaseSpaceShiftX",po::value<decltype(meshshiftx)>(
             &meshshiftx)->default_value(0),
            "Shift grid by X mesh points")
        ("PhaseSpaceShiftY",po::value<decltype(meshshifty)>(
             &meshshifty)->default_value(0),
            "Shift grid by Y mesh points")
        ("RenormalizeCharge",po::value<decltype(renormalize)>(
             &renormalize)->default_value(0),
            ">0: renormalize charge every n-th simulation step\n"
            " 0: do just one initial renormalization\n"
            "<0: no renormalization")
        ("FPType", po::value<decltype(fptype)>(
             &fptype)->default_value(3),
            "Used implementation for Fokker-Planck term\n"
            " 0: No Fokker-Planck term\n"
            " 1: Only damping\n"
            " 2: Only diffusion\n"
            " 3: Full")
        ("FPTrack", po::value<decltype(fptrack)>(
             &fptrack)->default_value(3),
            "Used implementation for FP term in tracking\n"
            " 0: No Fokker-Planck term\n"
            " 1: Approximation overweighting damping\n"
            " 2: Approximation overweighting diffusion\n"
            " 3: Stochastic")
        ("GridSize,s", po::value<decltype(meshsize)>(
             &meshsize)->default_value(256),
            "Number of mesh points per dimension")
        ("rotations,T", po::value<decltype(rotations)>(
             &rotations)->default_value(5),
            "Simulated time (in number of synchrotron periods)")
        ("derivation",po::value<decltype(deriv_type)>(
             &deriv_type)->default_value(4u),
            "Number of grid points to be used to numerically find derivative")
        ("InterpolationPoints",po::value<decltype(interpol_type)>(
             &interpol_type)->default_value(4u),
            "Number of grid points to be used for interpolation")
        ("InterpolateClamped",po::value<decltype(interpol_clamp)>(
             &interpol_clamp)->default_value(false),
            "Restrict result of interpolation to the"
            "values of the neighboring grid points")
    ;
    _compatopts_ignore.add_options()
        ("HaissinskiIterations",po::value<decltype(_hi)>(
             &_hi)->default_value(0),
            "(currently ignored)")
        ("InitialDistParam",po::value<decltype(_hi)>(
             &_hi)->default_value(0),
            "(currently ignored)")
        ("RotationType", po::value<decltype(rotationtype)>(
             &rotationtype),
            "compatibility (ignored)")
        ("SaveSourceMap", po::value<decltype(_savesourcemap)>(
             &_savesourcemap),
            "compatibility (ignored)")
    ;
    _compatopts_alias.add_options()
        ("RFVoltage", po::value<decltype(V_RF)>(
             &V_RF),
            "compatibility naming for AcceleratingVoltage")
        ("SyncFreq", po::value<decltype(f_s)>(
             &f_s),
            "(compatibility naming for SynchrotronFrequency)")
        ("steps", po::value<decltype(steps_per_Ts)>(
             &steps_per_Ts),
            "(compatibility naming for StepsPerTs)")
    ;
    _compatopts.add(_compatopts_ignore);
    _compatopts.add(_compatopts_alias);
    _cfgfileopts.add(_physopts);
    _cfgfileopts.add(_programopts_file);
    _cfgfileopts.add(_simulopts);
    _cfgfileopts.add(_compatopts);
    _commandlineopts.add(_proginfoopts);
    _commandlineopts.add(_programopts_cli);
    _commandlineopts.add(_simulopts);
    _commandlineopts.add(_physopts);
    _visibleopts.add(_commandlineopts);
}

/**
 * @brief parse parse command line parameters
 * @param argc argument count (handed over from main)
 * @param argv argument values (handed over from main)
 * @return Run simulation?
 *
 * Function will return true, except in cases where it is undesired that the
 * simulation is run. This is the case when an error occurs or if one of the
 * following command line parameters is present:
 * - `help`
 * - `copyright`
 * - `version`
 * - `buildinfo`
 * - `cldev < 0`
 */
bool vfps::ProgramOptions::parse(int ac, char** av)
{
    po::store(po::parse_command_line(ac, av, _commandlineopts), _vm);
    po::notify(_vm);

    if (_vm.count("help")) {
        std::cout << _visibleopts << std::endl;
        return false;
    }
    if (_vm.count("copyright")) {
        std::cout << vfps::copyright_notice() << std::endl;
        return false;
    }
    if (_vm.count("version")) {
        std::cout << vfps::inovesa_version(false) << std::endl;
        return false;
    }
    if (_vm.count("buildinfo")) {
        std::cout << vfps::inovesa_version(true) << std::endl;
        return false;
    }
    #if INOVESA_USE_OPENCL == 1
    if (_cldevice < 0) {
        OCLH::listCLDevices();
        return false;
    }
    #endif // INOVESA_USE_OPENCL
    if (_configfile == "/dev/null") {
        _configfile.clear();
    } else if (!_configfile.empty()) {
        if (boost::filesystem::exists(_configfile) &&
            boost::filesystem::is_regular_file(_configfile) ) {
            std::ifstream ifs(_configfile.c_str());
            if (!ifs) {
                std::cout << "Cannot open config file: " << _configfile
                          << std::endl;
                return false;
            } else {
                std::string message = "Loading configuration from \""
                                     + _configfile + "\".";
                Display::printText(message);
                store(parse_config_file(ifs, _cfgfileopts), _vm);
                notify(_vm);
                if(_vm.count("SyncFreq")) {
                    _vm.at("SynchrotronFrequency").value()
                            = _vm["SyncFreq"].value();
                }
                notify(_vm);
            }
        } else if (_configfile != "default.cfg") {
            std::cout << "Config file \"" << _configfile
                      << "\" does not exist."<< std::endl;
            return false;
        }
    }
    if (_outfile == "/dev/null") {
        _outfile.clear();
    }
    if (_startdistfile == "/dev/null") {
        _startdistfile.clear();
    }
    #if INOVESA_USE_OPENCL == 0
    if (_cldevice > 0 && _vm.count("cldev")) {
        std::cout   << "Warning: Defined device for OpenCL "
                    << "but running Inovesa without OpenCL support."
                    << std::endl;
    }
    #endif

    return true;
}

/**
 * @brief save settings to a configutation file
 * @param fname file name to use for saving
 *
 * @todo ignore members of _compatopts
 */
void vfps::ProgramOptions::save(std::string fname)
{
    std::ofstream ofs(fname.c_str());

    ofs << "# " << vfps::inovesa_version() << std::endl;

    for ( auto it=_vm.begin(); it != _vm.end(); it++ ) {
        // currently, the _compatopts are ignored manually
        if(it->first == "HaissinskiIterations"
        || it->first == "InitialDistParam"
        || it->first == "RotationType"
        || it->first == "SyncFreq"
        || it->first == "steps"
        || it->first == "RFVoltage"
        || it->first == "run_anyway"
        || it->first == "SaveSourceMap"
        ){
            continue;
        } else
        if (it->first == "alpha0" && std::fpclassify(f_s) == FP_ZERO) {
            ofs << "alpha0=0" << std::endl;
            continue;
        } else
        if (!it->second.value().empty()) {
            if (it->second.value().type() == typeid(float)) {
                ofs << it->first << '='
                    << _vm[it->first].as<float>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(double)) {
                ofs << it->first << '='
                    << _vm[it->first].as<double>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(int32_t)) {
                ofs << it->first << '='
                    << _vm[it->first].as<int32_t>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(uint32_t)) {
                ofs << it->first << '='
                    << _vm[it->first].as<uint32_t>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(int64_t)) {
                ofs << it->first << '='
                    << _vm[it->first].as<int64_t>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(bool)) {
                ofs << it->first << '='
                    << _vm[it->first].as<bool>()
                    << std::endl;
            } else {
                std::string val;
                try {
                    if (it->first == "config") {
                        ofs << '#';
                    }
                    val = _vm[it->first].as<std::string>();
                    ofs << it->first << '='
                        << val
                        << std::endl;
                } catch(const boost::bad_any_cast &){}
            }
        }
    }
}

#if INOVESA_USE_HDF5 == 1
void vfps::ProgramOptions::save(vfps::HDF5File* file)
{
    for (po::variables_map::iterator it=_vm.begin(); it != _vm.end(); ++it ) {
        if (!it->second.value().empty()) {
            if (it->second.value().type() == typeid(float)) {
                float data = _vm[it->first].as<float>();
                file->addParameterToGroup("/Info/Parameters",it->first,
                                          H5::PredType::IEEE_F32LE,&data);
            } else if (it->second.value().type() == typeid(double)) {
                double data = _vm[it->first].as<double>();
                file->addParameterToGroup("/Info/Parameters",it->first,
                                          H5::PredType::IEEE_F64LE,&data);
            } else if (it->second.value().type() == typeid(uint32_t)) {
                uint32_t data = _vm[it->first].as<uint32_t>();
                file->addParameterToGroup("/Info/Parameters",it->first,
                                          H5::PredType::STD_U32LE,&data);
            } else if (it->second.value().type() == typeid(uint64_t)) {
                uint64_t data = _vm[it->first].as<uint64_t>();
                file->addParameterToGroup("/Info/Parameters",it->first,
                                          H5::PredType::STD_U32LE,&data);
            } else if (it->second.value().type() == typeid(int32_t)) {
                int32_t data = _vm[it->first].as<int32_t>();
                file->addParameterToGroup("/Info/Parameters",it->first,
                                          H5::PredType::STD_I32LE,&data);
            }
        }
    }
}
#endif // INOVESA_USE_HDF5
