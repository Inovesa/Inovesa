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

#include "IO/ProgramOptions.hpp"

vfps::ProgramOptions::ProgramOptions() :
    _configfile("default.cfg"),
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
    ;
    _physopts.add_options()
        ("alpha0", po::value<double>(&alpha0)->default_value(4e-3,"4e-3"),
            "Linear Momentum compaction factor (1)")
        ("alpha1", po::value<double>(&alpha1)->default_value(0),
            "Quadratic Momentum compaction factor (1)")
        ("alpha2", po::value<double>(&alpha2)->default_value(0),
            "Cubic Momentum compaction factor (1)")
        ("RevolutionFrequency,F",po::value<double>(&f0)->default_value(9e6,"9e6"),
            "Revolution frequency (Hz)")
        ("DampingTime,d", po::value<double>(&t_d)->default_value(0.001),
            "Damping time (s)")
        ("HarmonicNumber,H", po::value<double>(&H)->default_value(50),
            "Harmonic Number (1)")
        ("InitialDistFile,i", po::value<std::string>(&_startdistfile),
            "might be:\n"
             #ifdef INOVESA_USE_HDF5
             "\tInovesa result file (HDF5, .h5)\n"
             #endif // INOVESA_USE_HDF5
             #ifdef INOVESA_USE_PNG
             "\tgrayscale png (.png) file\n"
             #endif // INOVESA_USE_PNG
             "\ttext file (.txt) w/ particle coordinates")
        ("InitialDistParam",po::value<double>(&Fk)->default_value(0),
            "Parameter F(k) for generation of initial distribution")
        ("InitialDistZoom",po::value<double>(&zoom)->default_value(1),
            "Magnification for generation of initial distribution")
        ("BunchCurrent,I", po::value<double>(&I_b)->default_value(3e-3,"3e-3"),
            "Ring Current due to a single bunch (A)")
        ("BendingRadius,R", po::value<double>(&r_bend)->default_value(-1),
            "Bending radius of accelerator (m)\n"
            "negative: calculate from RevolutionFrequency")
        ("BeamEnergy,E", po::value<double>(&E_0)->default_value(1.3e9,"1.3e9"),
            "Beam energy (eV)")
        ("BeamEnergySpread,e", po::value<double>(&s_E)->default_value(4.7e-4,"4.7e-4"),
            "Natural energy spread (relative)")
        ("Impedance,Z", po::value<std::string>(&_impedancefile),
            "File containing impedance information.")
        ("VacuumGap,G", po::value<double>(&g)->default_value(0.03,"0.03"),
            "Full height of vacuum chamber (m)\n"
            "<0: free space CSR\n"
            " 0: no CSR\n"
            ">0: parallel plates CSR")
        ("CutoffFreq", po::value<double>(&f_c)->default_value(23e9,"23e9"),
            "Beamline cutoff frequency (Hz)")
        ("AcceleratingVoltage,V",
            po::value<double>(&V_RF)->default_value(1e6,"1e6"),
            "Accelerating Voltage (V) for one revolution")
        ("WakeFunction,w", po::value<std::string>(&_wakefile),
            "File containing wake function.")
    ;
    _programopts_file.add_options()
        ("cldev", po::value<int>(&_cldevice)->default_value(1),
            "OpenCL device to use\n('-1' lists available devices)")
        ("ForceOpenGLVersion", po::value<int>(&_glversion)->default_value(2),
            "Force OpenGL version")
        ("gui,g", po::value<bool>(&_showphasespace)->default_value(true),
            "Show phase space view")
        ("output,o",
            po::value<std::string>(&_outfile),
            "name of file to safe results.")
        ("SavePhaseSpace",
            po::value<bool>(&_savephasespace)->default_value(false),
            "save every outstep's phase space to HDF5 file")
        ("tracking",
            po::value<std::string>(&_trackingfile)->default_value(""),
            "file containing starting positions (grid points)"
            "of particles to be (pseudo-) tracked")
        ("verbose,v", po::value<bool>(&_verbose)->default_value(false),
            "print information more detailed")
    ;
    _programopts_cli.add_options()
        #ifdef INOVESA_USE_CL
        ("cldev", po::value<int>(&_cldevice)->default_value(0),
            "OpenCL device to use\n('-1' lists available devices)")
        #endif // INOVESA_USE_CL
        ("config,c", po::value<std::string>(&_configfile),
            "name of a file containing a configuration.")
        ("ForceOpenGLVersion", po::value<int>(&_glversion)->default_value(2),
            "Force OpenGL version")
        ("gui,g", po::value<bool>(&_showphasespace)->default_value(true),
            "Show phase space view")
        ("output,o",
            po::value<std::string>(&_outfile),
            "name of file to safe results.")
        ("SavePhaseSpace",
            po::value<bool>(&_savephasespace)->default_value(false),
            "save every outstep's phase space to HDF5 file")
        ("tracking",
            po::value<std::string>(&_trackingfile)->default_value(""),
            "file containing starting positions (grid points)"
            "of particles to be (pseudo-) tracked")
        ("verbose,v", "print information more detailed")
    ;
    _simulopts.add_options()
        ("steps,N", po::value<uint32_t>(&steps)->default_value(1000),
            "Steps for one synchrotron period")
        ("outstep,n", po::value<uint32_t>(&outsteps)->default_value(100),
            "Save results every n steps.")
        ("padding,p", po::value<double>(&padding)->default_value(8.0),
            "Factor for zero padding of bunch profile")
        ("PhaseSpaceSize,P", po::value<double>(&pq_size)->default_value(12),
            "Size of phase space")
        ("PhaseSpaceShiftX",po::value<double>(&meshshiftx)->default_value(0),
            "Shift grid by X mesh points")
        ("PhaseSpaceShiftY",po::value<double>(&meshshifty)->default_value(0),
            "Shift grid by Y mesh points")
        ("RenormalizeCharge",po::value<bool>(&renormalize)->default_value(false),
            "Reset charge every outstep")
        ("RotationType", po::value<uint32_t>(&rotationtype)->default_value(2),
            "Used implementation for rotation\n"
            " 0: Standard rotation without source map\n"
            " 1: Standard rotation with source map\n"
            " 2: Manhattan rotation")
        ("GridSize,s", po::value<uint32_t>(&meshsize)->default_value(256),
            "Number of mesh points per dimension")
        ("rotations,T", po::value<double>(&rotations)->default_value(5),
            "Simulated time (in number of synchrotron periods)")
        ("derivation",po::value<uint32_t>(&deriv_type)->default_value(4u),
            "Number of grid points to be used to numerically find derivative")
        ("InterpolationPoints",po::value<uint32_t>(&interpol_type)->default_value(4u),
            "Number of grid points to be used for interpolation")
        ("InterpolateClamped",po::value<bool>(&interpol_clamp)->default_value(false),
            "Restrict result of interpolation to the values of the neighboring grid points")
    ;
    _cfgfileopts.add(_physopts);
    _cfgfileopts.add(_programopts_file);
    _cfgfileopts.add(_simulopts);
    _commandlineopts.add(_proginfoopts);
    _commandlineopts.add(_programopts_cli);
    _commandlineopts.add(_simulopts);
    _commandlineopts.add(_physopts);
    _visibleopts.add(_commandlineopts);


    std::stringstream timestamp;
    timestamp << time(nullptr);
    _outfile = "inovesa-result_" + timestamp.str() + ".h5";
}

bool vfps::ProgramOptions::parse(int ac, char** av)
{
    po::store(po::parse_command_line(ac, av, _commandlineopts), _vm);
    po::notify(_vm);

    if (_vm.count("verbose")) {
        _verbose = true;
    }
    if (_vm.count("help")) {
        std::cout << _visibleopts << std::endl;
        return false;
    }
    if (_vm.count("copyright")) {
        std::cout << vfps::copyright_notice() << std::endl;
        return false;
    }
    if (_vm.count("version")) {
        std::cout << vfps::inovesa_version() << std::endl;
        return false;
    }
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
        }
    } else if (_configfile != "default.cfg") {
        std::cout << "Config file \"" << _configfile
                  << "\" does not exist."<< std::endl;
    }
    #ifndef INOVESA_USE_CL
    if (_vm.count("cldev")) {
        std::cout    << "Warning: Defined device for OpenCL "
                    << "but running Inovesa without OpenCL support."
                    << std::endl;
    }
    #endif

    return true;
}

void vfps::ProgramOptions::save(std::string fname)
{
    std::ofstream ofs(fname.c_str());

    ofs << "# " << vfps::inovesa_version() << std::endl;

    for (po::variables_map::iterator it=_vm.begin(); it != _vm.end(); it++ ) {
        if (!it->second.value().empty()) {
            if (it->second.value().type() == typeid(double)) {
                ofs << it->first << '='
                    << _vm[it->first].as<double>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(uint32_t)) {
                ofs << it->first << '='
                    << _vm[it->first].as<uint32_t>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(int32_t)) {
                ofs << it->first << '='
                    << _vm[it->first].as<int32_t>()
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

void vfps::ProgramOptions::save(vfps::HDF5File* file)
{
    for (po::variables_map::iterator it=_vm.begin(); it != _vm.end(); it++ ) {
        if (!it->second.value().empty()) {
            if (it->second.value().type() == typeid(double)) {
                double data = _vm[it->first].as<double>();
                file->addParameterToGroup("/Info/Parameters",it->first,
                                          H5::PredType::IEEE_F64LE,&data);
            } else if (it->second.value().type() == typeid(uint32_t)) {
                uint32_t data = _vm[it->first].as<uint32_t>();
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
