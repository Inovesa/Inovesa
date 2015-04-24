#ifndef PROGRAMOPTIONS_HPP
#define PROGRAMOPTIONS_HPP

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <string>

#include "defines.hpp"

namespace po = boost::program_options;

namespace vfps
{

class ProgramOptions
{
public:
	ProgramOptions();

	bool parse(int argc, char** argv);

	inline unsigned int getCLDevice() const
		{ return _cldevice; }

	inline std::string getStartDistPNG() const
		{ return _startdistpng; }

private:
	unsigned int _cldevice;

	std::string _startdistpng;

	std::string _configfile;

private:
	po::options_description _commandlineopts;

	po::options_description _configopts;

	po::options_description _cfgfileopts;

	po::options_description _genericopts;

	po::options_description _hiddenopts;

	po::options_description _visibleopts;

	po::variables_map _vm;
};

}

#endif // PROGRAMOPTIONS_HPP
