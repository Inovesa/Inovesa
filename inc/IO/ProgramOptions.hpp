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

	inline std::string getOutFile() const
		{ return _outfile; }

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

	std::string _outfile;

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
