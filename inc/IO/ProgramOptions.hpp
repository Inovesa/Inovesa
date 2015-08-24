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
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "defines.hpp"

namespace po = boost::program_options;

namespace vfps
{

class ProgramOptions
{
public:
	ProgramOptions();

	bool parse(int argc, char** argv);

	inline int getCLDevice() const
		{ return _cldevice; }

	inline std::string getImpedanceFile() const
		{ return _impedancefile; }

	inline std::string getOutFile() const
		{ return _outfile; }

	inline bool showPhaseSpace() const
		{ return _showphasespace; }

	inline std::string getStartDistFile() const
		{ return _startdistfile; }

public:
	inline unsigned int getOutSteps() const
		{ return outsteps; }

	inline unsigned int getPadding() const
		{ return padding; }

	inline unsigned int getSteps() const
		{ return steps; }

	inline float getNRotations() const
		{ return rotations; }


public:
	inline double getBendingRadius() const
		{ return r_bend; }

	inline double getNaturalBunchLength() const
		{ return s_0; }

	inline double getSyncFreq() const
		{ return f_s; }

	inline double getDampingTime() const
		{ return t_d; }

	inline double getPhaseSpaceSize() const
		{ return pq_max; }

private: // program parameters
	int _cldevice;

	std::string _impedancefile;

	std::string _outfile;

	bool _showphasespace;

	std::string _startdistfile;

	std::string _configfile;

private: // simulation parameters
	unsigned int outsteps;
	unsigned int padding;
	unsigned int steps;
	float rotations;

private: // phsical parameters
	double f_s;
	double t_d;
	double pq_max;
	double r_bend;
	double s_0;

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
