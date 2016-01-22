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

#include "IO/ProgramOptions.hpp"

vfps::ProgramOptions::ProgramOptions() :
    _cldevice(1),
    _startdistfile(""),
    _configfile("default.cfg"),
    _wakefile(""),
    meshsize(256),
    outsteps(100),
    padding(0),
    pq_max(5.0),
    steps(4000),
    rotations(1),
    rotmapsize(-1),
    E_0(1.3e9),
    Fk(0),
    f_c(23e9),
    f_s(8.5e3),
    f_rev(2.7e6),
    h(1),
    I_b(1),
    t_d(0.01),
    r_bend(1.0),
    s_E(4.7e-4),
    V_RF(1e6),
    _physopts("Physical Parameters for Simulation"),
    _proginfoopts("Program Information"),
    _programopts_cli("General Program Parameters"),
    _simulopts("Non-Physical Parameters for Simulation"),
    _visibleopts("Possible Parameters")
{
    _proginfoopts.add_options()
        ("help,h", "print help message")
        ("version,v", "print version string")
    ;
    _physopts.add_options()
        ("RevolutionFrequency,F", po::value<double>(&f_rev),
            "Revolution frequency (Hz)")
        ("SyncFreq,f", po::value<double>(&f_s),"Synchrotron frequency (Hz)")
        ("DampingTime,d", po::value<double>(&t_d),"Damping time (s)")
        ("HarmonicNumber,H", po::value<double>(&h),
            "Harmonic Number (1)")
        ("InitialDistFile,D", po::value<std::string>(&_startdistfile),
            "might be:\n"
            "\tgrayscale png (.png) file containing initial particle density\n"
            "\ttext file (.txt) containing normalized z/delta of particles")
        ("InitialDistParam,K",po::value<double>(&Fk),
            "Parameter F(k) of initial distribution")
        ("BunchCurrent,I", po::value<double>(&I_b),
            "Ring Current due to a single electron bunch (A)")
        ("BendingRadius,R", po::value<double>(&r_bend),
            "Bending radius of accelerator (m)")
        ("BeamEnergy,E", po::value<double>(&E_0),
            "Beam energy (GeV)")
        ("BeamEnergySpread,e", po::value<double>(&s_E),
            "Natural energy spread (relative)")
        ("Impedance,Z", po::value<std::string>(&_impedancefile),
            "File containing impedance information. Might be:\n"
            "\ttext file (.txt) containing wavenumber Re(Z) Im(Z)")
        ("CutoffFreq", po::value<double>(&f_c),"Beamline cutoff frequency (Hz)")
        ("RFVoltage,V", po::value<double>(&V_RF),
            "Accelerating Voltage (V)")
        ("WakeFunction,w", po::value<std::string>(&_wakefile),
            "File containing information on wake function.")
    ;
    _programopts_file.add_options()
        ("cldev", po::value<int>(&_cldevice)->default_value(1),
            "OpenCL device to use ('0' lists available devices)")
        ("gui,g", po::value<bool>(&_showphasespace)->default_value(true),
            "Show phase space view")
        ("output,o",
            po::value<std::string>(&_outfile),
            "name of file to safe results.")
    ;
    _programopts_cli.add_options()
        #ifdef INOVESA_USE_CL
        ("cldev", po::value<int>(&_cldevice)->default_value(1),
            "OpenCL device to use ('0' lists available devices)")
        #endif // INOVESA_USE_CL
        ("config,c", po::value<std::string>(&_configfile),
            "name of a file containing a configuration.")
        ("gui,g", po::value<bool>(&_showphasespace)->default_value(true),
            "Show phase space view")
        ("output,o",
            po::value<std::string>(&_outfile),
            "name of file to safe results.")
    ;
    _simulopts.add_options()
        ("steps,N", po::value<unsigned int>(&steps),
            "Number of steps for one synchrotron period (delta t=1/(N*f_s))")
        ("outstep,n", po::value<unsigned int>(&outsteps),
            "Save results/ update phase space view every n steps")
        ("padding,p", po::value<unsigned int>(&padding),
            "Factor to use for zero padding of bunch profile, 0/1: no padding")
        ("PhaseSpaceSize,P", po::value<double>(&pq_max),
            "Size of phase space (maximum sigma_{p/q})")
        ("ForceRotMapSize", po::value<int>(&rotmapsize)->default_value(-1),
            "Size of rotation map (-1=guess, 0=none, 1=fill, 2=half)")
        ("MeshSize,s", po::value<unsigned int>(&meshsize),
            "Size of phase space mesh (number of p/q mesh points)")
        ("rotations,T", po::value<double>(&rotations),
            "Simulated time (in number of synchrotron periods)")
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
    _outfile = "result_" + timestamp.str() + ".h5";
}

bool vfps::ProgramOptions::parse(int ac, char** av)
{
    po::store(po::parse_command_line(ac, av, _commandlineopts), _vm);
    po::notify(_vm);

    if (_vm.count("help")) {
        std::cout << _visibleopts << std::endl;
        return false;
    }
    if (_vm.count("version")) {
        std::cout << "Inovesa v"
                  << INOVESA_VERSION_RELEASE << '.'
                  << INOVESA_VERSION_MINOR << '.'
                  << INOVESA_VERSION_FIX;
        if (std::string(GIT_BRANCH) != "stable") {
            std::cout << " (Branch: " GIT_BRANCH ")";
        }
        std::cout << std::endl;
        return false;
    }
    std::ifstream ifs(_configfile.c_str());
    if (!ifs) {
        std::cout << "Cannot open config file: " << _configfile
                  << std::endl;
        return false;
    } else {
        store(parse_config_file(ifs, _cfgfileopts), _vm);
        notify(_vm);
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

    ofs << "#Inovesa v"
        << INOVESA_VERSION_RELEASE << '.'
        << INOVESA_VERSION_MINOR << '.'
        << INOVESA_VERSION_FIX;
    if (std::string(GIT_BRANCH) != "stable") {
        ofs << " (Branch: " GIT_BRANCH ")";
    }
    ofs << std::endl;

    for (po::variables_map::iterator it=_vm.begin(); it != _vm.end(); it++ ) {
        if (!it->second.value().empty() && !_vm[it->first].defaulted()) {
            if (it->second.value().type() == typeid(double)) {
                ofs << it->first << '='
                    << _vm[it->first].as<double>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(unsigned int)) {
                ofs << it->first << '='
                    << _vm[it->first].as<unsigned int>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(int)) {
                ofs << it->first << '='
                    << _vm[it->first].as<int>()
                    << std::endl;
            } else if (it->second.value().type() == typeid(bool)) {
                ofs << it->first << '='
                    << _vm[it->first].as<bool>()
                    << std::endl;
            } else {
                std::string val;
                try {
                    val = _vm[it->first].as<std::string>();
                    ofs << it->first << '='
                        << val
                        << std::endl;
                } catch(const boost::bad_any_cast &){}
            }
        }
    }
}
