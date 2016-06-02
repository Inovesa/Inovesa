/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Algorithms   *
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

#ifndef PROGRAMOPTIONS_HPP
#define PROGRAMOPTIONS_HPP

#include <boost/program_options.hpp>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "defines.hpp"
#include "SM/FokkerPlanckMap.hpp"

namespace po = boost::program_options;

namespace vfps
{

class ProgramOptions
{
public:
    ProgramOptions();

    bool parse(int argc, char** argv);

    void save(std::string fname);

public:
    inline int getCLDevice() const
        { return _cldevice; }

    inline std::string getImpedanceFile() const
        { return _impedancefile; }

    inline std::string getOutFile() const
        { return _outfile; }

    inline int getOpenGLVersion() const
        { return _glversion; }

    inline bool showPhaseSpace() const
        { return _showphasespace; }

    inline std::string getStartDistFile() const
        { return _startdistfile; }

    inline bool getVerbosity() const
        { return _verbose; }

    inline std::string getWakeFile() const
        { return _wakefile; }

public:
    inline unsigned int getMeshSize() const
        { return meshsize; }

    inline unsigned int getOutSteps() const
        { return outsteps; }

    inline unsigned int getPadding() const
        { return padding; }

    inline unsigned int getSteps() const
        { return steps; }

    inline float getNRotations() const
        { return rotations; }

    inline double getPhaseSpaceSize() const
        { return pq_size; }

    inline double getPSShiftX() const
        { return meshshiftx; }

    inline double getPSShiftY() const
        { return meshshifty; }

    inline int getRotationType() const
        { return rotationtype; }

    inline unsigned int getDerivationType() const
        { return deriv_type; }

    inline unsigned int getInterpolationPoints() const
        { return interpol_type; }

    inline bool getInterpolationBound() const
        { return interpol_clamp; }

public:
    inline double getBeamEnergy() const
        { return E_0; }

    inline double getBendingRadius() const
        { return r_bend; }

    inline double getBunchCurrent() const
        { return I_b; }

    inline double getCutoffFrequency() const
        { return f_c; }

    inline double getEnergySpread() const
        { return s_E; }

    inline double getHarmonicNumber() const
        { return H; }

    inline double getRevolutionFrequency() const
        { return f_rev; }

    inline double getRFVoltage() const
        { return V_RF; }

    inline double getStartDistParam() const
        { return Fk; }

    inline double getSyncFreq() const
        { return f_s; }

    inline double getDampingTime() const
        { return t_d; }

    inline double getVacuumChamberHeight() const
        { return h; }

private: // program parameters
    int _cldevice;

    std::string _impedancefile;

    std::string _outfile;

    bool _showphasespace;

    std::string _startdistfile;

    std::string _configfile;

    int _glversion;

    bool _verbose;

    std::string _wakefile;

private: // simulation parameters
    unsigned int meshsize;
    unsigned int outsteps;
    unsigned int padding;
    double pq_size;
    double meshshiftx;
    double meshshifty;
    unsigned int steps;
    double rotations;
    unsigned int rotationtype;
    unsigned int deriv_type;
    unsigned int interpol_type;
    bool interpol_clamp;

private: // phsical parameters
    double E_0;
    double Fk;
    double f_c;
    double f_s;
    double f_rev;
    double h;
    double H;
    double I_b;
    double t_d;
    double r_bend;
    double s_E;
    double V_RF;

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
