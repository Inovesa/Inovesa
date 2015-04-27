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

public:
	inline unsigned int getOutSteps() const
		{ return outsteps; }

	inline unsigned int getSteps() const
		{ return steps; }

	inline float getNRotations() const
		{ return rotations; }


public:
	inline double getSyncFreq()
		{ return f_s; }

	inline double getDampingTime()
		{ return t_d; }

private: // program parameters
	unsigned int _cldevice;

	std::string _startdistpng;

	std::string _configfile;

private: // simulation parameters
	unsigned int outsteps;
	unsigned int steps;
	float rotations;

private: // phsical parameters
	double f_s;
	double t_d;

private:
	po::options_description _cfgfileopts;

	po::options_description _commandlineopts;

	po::options_description _hiddenopts;

	po::options_description _physopts;

	po::options_description _proginfoopts;

	po::options_description _programopts_cli;

	po::options_description _programopts_file;

	po::options_description _simulopts;

	po::options_description _visibleopts;

	po::variables_map _vm;
};

}

#endif // PROGRAMOPTIONS_HPP
