/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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
#include "MessageStrings.hpp"
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

    /**
     * @brief save
     * @param fname
     *
     * @todo ignore members of _compatopts
     */
    void save(std::string fname);

    #ifdef INOVESA_USE_HDF5
    void save(HDF5File* file);
    #endif // INOVESA_USE_HDF5

public:
    inline auto getCLDevice() const
        { return _cldevice; }

    inline auto getImpedanceFile() const
        { return _impedancefile; }

    inline auto getOutFile() const
        { return _outfile; }

    inline auto getSavePhaseSpace() const
        { return _savephasespace; }

    inline auto getSaveSourceMap() const
        { return _savesourcemap; }

    #ifdef INOVESA_USE_OPENGL
    inline auto getOpenGLVersion() const
        { return _glversion; }

    inline auto showPhaseSpace() const
        { return _showphasespace; }
    #endif // INOVESA_USE_OPENGL

    inline auto getStartDistFile() const
        { return _startdistfile; }

    inline auto getStartDistStep() const
        { return _startdiststep; }

    inline auto getParticleTracking() const
        { return _trackingfile; }

    inline auto getVerbosity() const
        { return _verbose; }

    inline auto getForceRun() const
        { return _forcerun; }

    inline auto getWakeFile() const
        { return _wakefile; }

public:
    inline auto getGridSize() const
        { return meshsize; }

    inline auto getOutSteps() const
        { return outsteps; }

    inline auto getPadding() const
        { return padding; }

    inline auto getRoundPadding() const
        { return roundpadding; }

    inline auto getStepsPerTsync() const
        { return steps_per_Ts; }

    inline auto getStepsPerTrev() const
        { return steps_per_Trev; }

    inline auto getNRotations() const
        { return rotations; }

    inline auto getPhaseSpaceSize() const
        { return pq_size; }

    inline auto getPSShiftX() const
        { return meshshiftx; }

    inline auto getPSShiftY() const
        { return meshshifty; }

    inline auto getRenormalizeCharge() const
        { return renormalize; }

    inline auto getFPTrack() const
        { return fptrack; }

    inline auto getFPType() const
        { return fptype; }

    inline auto getDerivationType() const
        { return deriv_type; }

    inline auto getInterpolationPoints() const
        { return interpol_type; }

    inline auto getInterpolationClamped() const
        { return interpol_clamp; }

public:
    inline auto getAlpha0() const
        { return alpha0; }

    inline auto getAlpha1() const
        { return alpha1; }

    inline auto getAlpha2() const
        { return alpha2; }

    inline auto getRFAmplitudeSpread() const
        { return rf_amplitude_spread; }

    inline auto getRFPhaseSpread() const
        { return rf_phase_spread; }

    inline auto getRFPhaseModAmplitude() const
        { return rf_phase_mod_amplitude; }

    inline auto getRFPhaseModFrequency() const
        { return rf_phase_mod_frequency; }

    inline auto getBeamEnergy() const
        { return E_0; }

    inline auto getBendingRadius() const
        { return r_bend; }

    inline auto getBunchCurrent() const
        { return I_b; }

    inline auto getCutoffFrequency() const
        { return f_c; }

    inline auto getEnergySpread() const
        { return s_E; }

    inline auto getHaissinskiIterations() const
        { return _hi; }

    inline auto getHarmonicNumber() const
        { return H; }

    inline auto getRevolutionFrequency() const
        { return f0; }

    inline auto getRFVoltage() const
        { return V_RF; }

    inline auto getStartDistZoom() const
        { return zoom; }

    inline auto getSyncFreq() const
        { return f_s; }

    inline auto getDampingTime() const
        { return t_d; }

    inline auto getVacuumChamberGap() const
        { return g; }

    inline auto getUseCSR() const
        { return use_csr; }

    inline auto getLinearRF() const
        { return linearRF; }

    inline auto getCollimatorRadius() const
        { return collimator; }

    inline auto getWallConductivity() const
        { return s_c; }

    inline auto getWallSusceptibility() const
        { return xi_wall; }

private: // program parameters
    int32_t _cldevice;

    std::string _impedancefile;

    std::string _outfile;

    uint32_t _savephasespace;

    bool _savesourcemap;

    bool _showphasespace;

    std::string _startdistfile;

    std::string _trackingfile;

    int64_t _startdiststep;

    std::string _configfile;

    int32_t _glversion;

    bool _verbose;

    bool _forcerun;

    std::string _wakefile;

private: // simulation parameters
    uint32_t meshsize;
    uint32_t outsteps;
    double padding;
    bool roundpadding;
    double pq_size;
    double meshshiftx;
    double meshshifty;
    uint32_t steps_per_Ts;
    double steps_per_Trev;
    int32_t renormalize;
    double rotations;
    uint32_t fptype;
    uint32_t fptrack;
    uint32_t rotationtype;
    uint32_t deriv_type;
    uint32_t interpol_type;
    bool interpol_clamp;

private: // phsical parameters
    double alpha0;
    double alpha1;
    double alpha2;

    double rf_phase_spread;
    double rf_amplitude_spread;

    double rf_phase_mod_amplitude;
    double rf_phase_mod_frequency;

    double E_0;

    uint32_t _hi;

    /**
     * @brief zoom initial distribution
     */
    double zoom;
    double f_c;
    double f_s;
    double f0;
    double g;
    double collimator;
    double s_c;
    double xi_wall;
    double H;
    double I_b;
    double t_d;
    double r_bend;
    double s_E;
    double V_RF;
    bool linearRF;

    bool use_csr;

private:
    po::options_description _compatopts;

    po::options_description _cfgfileopts;

    po::options_description _commandlineopts;

    po::options_description _hiddenopts;

    po::options_description _physopts;

    po::options_description _proginfoopts;

    po::options_description _programopts_cli;

    po::options_description _programopts_file;

    /**
     * @brief _compatopts options from old config files
     */
    po::options_description _compatopts_alias;

    po::options_description _compatopts_ignore;

    po::options_description _simulopts;

    po::options_description _visibleopts;

    po::variables_map _vm;
};

}

#endif // PROGRAMOPTIONS_HPP
