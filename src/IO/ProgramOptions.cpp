#include "IO/ProgramOptions.hpp"

vfps::ProgramOptions::ProgramOptions() :
	_cldevice(0),
	_startdistpng("start.png"),
	_configopts("Configuration"),
	_genericopts("Generic options"),
	_visibleopts("Possible Options")
{
	_genericopts.add_options()
		("help,h", "print help message")
		("version,v", "print version string")
	;
	_configopts.add_options()
		#ifdef INOVESA_USE_CL
		("cldev", po::value<unsigned int>(&_cldevice)->default_value(0),
			"OpenCL device to use ('0' lists available devices)")
		#endif // INOVESA_USE_CL
		("start-dist,S", po::value<std::string>(&_startdistpng),
			"grayscale png file containing initial particle density")
		("config,c", po::value<std::string>(&_configfile),
			"name of a file containing a configuration.")
	;
	_commandlineopts.add(_genericopts).add(_configopts);
	_cfgfileopts.add(_configopts);
	_visibleopts.add(_genericopts).add(_configopts);
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
	if (_vm.count("config")) {
		std::ifstream ifs(_configfile.c_str());
		if (!ifs) {
			std::cout << "Cannot open config file: " << _configfile << "\n";
			return false;
		} else {
			store(parse_config_file(ifs, _configopts), _vm);
			notify(_vm);
		}
	}

	return true;
}
