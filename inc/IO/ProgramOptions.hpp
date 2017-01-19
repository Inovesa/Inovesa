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

#ifndef PROGRAMOPTIONS_HPP
#define PROGRAMOPTIONS_HPP

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "defines.hpp"
#include "SM/FokkerPlanckMap.hpp"
#ifdef INOVESA_USE_HDF5
#include "IO/HDF5File.hpp"
#endif

namespace po = boost::program_options;

namespace vfps
{

class ProgramOptions
{
public:
    ProgramOptions();

    bool parse(int argc, char** argv);

    void save(std::string fname);

#ifdef INOVESA_USE_HDF5
    void save(HDF5File* file);
#endif // INOVESA_USE_HDF5

public:
    inline int getCLDevice() const
        { return _cldevice; }

    inline std::string getImpedanceFile() const
        { return _impedancefile; }

    inline std::string getOutFile() const
        { return _outfile; }

    inline bool getSavePhaseSpace() const
        { return _savephasespace; }

    inline int getOpenGLVersion() const
        { return _glversion; }

    inline bool showPhaseSpace() const
        { return _showphasespace; }

    inline std::string getStartDistFile() const
        { return _startdistfile; }

    inline int64_t getStartDistStep() const
        { return _startdiststep; }

    inline bool getVerbosity() const
        { return _verbose; }

    inline std::string getWakeFile() const
        { return _wakefile; }

public:
    inline uint32_t getMeshSize() const
        { return meshsize; }

    inline uint32_t getOutSteps() const
        { return outsteps; }

    inline double getPadding() const
        { return padding; }

    inline uint32_t getSteps() const
        { return steps; }

    inline float getNRotations() const
        { return rotations; }

    inline double getPhaseSpaceSize() const
        { return pq_size; }

    inline double getPSShiftX() const
        { return meshshiftx; }

    inline double getPSShiftY() const
        { return meshshifty; }

    inline uint32_t getRenormalizeCharge() const
        { return renormalize; }

    inline int getRotationType() const
        { return rotationtype; }

    inline uint32_t getDerivationType() const
        { return deriv_type; }

    inline uint32_t getInterpolationPoints() const
        { return interpol_type; }

    inline bool getInterpolationBound() const
        { return interpol_clamp; }

public:
    inline double getAlpha0() const
        { return alpha0; }

    inline double getAlpha1() const
        { return alpha1; }

    inline double getAlpha2() const
        { return alpha2; }

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
        { return f0; }

    inline double getRFVoltage() const
        { return V_RF; }

    inline double getStartDistParam() const
        { return Fk; }

    inline double getStartDistZoom() const
        { return zoom; }

    inline double getSyncFreq() const
        { return f_s; }

    inline double getDampingTime() const
        { return t_d; }

    inline double getVacuumChamberGap() const
        { return g; }

private: // program parameters
    int _cldevice;

    std::string _impedancefile;

    std::string _outfile;

    bool _savephasespace;

    bool _showphasespace;

    std::string _startdistfile;

    int64_t _startdiststep;

    std::string _configfile;

    int _glversion;

    bool _verbose;

    std::string _wakefile;

private: // simulation parameters
    uint32_t meshsize;
    uint32_t outsteps;
    double padding;
    double pq_size;
    double meshshiftx;
    double meshshifty;
    uint32_t steps;
    uint32_t renormalize;
    double rotations;
    uint32_t rotationtype;
    uint32_t deriv_type;
    uint32_t interpol_type;
    bool interpol_clamp;

private: // phsical parameters
    double alpha0;
    double alpha1;
    double alpha2;

    double E_0;
    double Fk;

    /**
     * @brief zoom initial distribution
     */
    double zoom;
    double f_c;
    double f_s;
    double f0;
    double g;
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
