/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2014-2015: Patrik Schönfeldt                                 */
/*                                                                            */
/* This file is part of Inovesa.                                              */
/* Inovesa is free software: you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* Inovesa is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           */
/******************************************************************************/

#include "IO/ProgramOptions.hpp"

vfps::ProgramOptions::ProgramOptions() :
	_cldevice(0),
	_startdistpng("start.png"),
	_configfile("default.cfg"),
	outsteps(100),
	steps(4000),
	rotations(1),
	f_s(8.5e3),
	t_d(0.01),
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
		("syncfreq,F", po::value<double>(&f_s),"Syncrotron frequency")
		("tdamp,T", po::value<double>(&t_d),"Damping time")
		("initial-dist,I", po::value<std::string>(&_startdistpng),
			"grayscale png file containing initial particle density")
	;
	_programopts_file.add_options()
		("cldev", po::value<unsigned int>(&_cldevice)->default_value(0),
			"OpenCL device to use ('0' lists available devices)")
		("gui", po::value<bool>(&_showphasespace)->default_value(true),
			"Show phase space view")
		("output,o",
			po::value<std::string>(&_outfile),
			"name of file to safe resuults.")
	;
	_programopts_cli.add_options()
		#ifdef INOVESA_USE_CL
		("cldev", po::value<unsigned int>(&_cldevice)->default_value(0),
			"OpenCL device to use ('0' lists available devices)")
		#endif // INOVESA_USE_CL
		("config,c", po::value<std::string>(&_configfile),
			"name of a file containing a configuration.")
		("gui,g", po::value<bool>(&_showphasespace)->default_value(true),
			"Show phase space view")
		("output,o",
			po::value<std::string>(&_outfile),
			"name of file to safe resuults.")
	;
	_simulopts.add_options()
		("steps,N", po::value<unsigned int>(&steps),
			"Number of steps for one synchrotron period (delta t=1/f_s)")
		("outstep,n", po::value<unsigned int>(&outsteps),
			"Save results every n steps")
		("rotations,R", po::value<float>(&rotations),
			"Number of totations to do")
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
				  << INOVESA_VERSION_FIX
				  << std::endl;
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
		std::cout	<< "Warning: Defined device for OpenCL "
					<< "but running Inovesa without OpenCL support."
					<< std::endl;
	}
	#endif

	return true;
}
