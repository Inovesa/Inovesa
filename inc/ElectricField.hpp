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

#ifndef ELECTRICFIELD_HPP
#define ELECTRICFIELD_HPP

#include <algorithm>
#include <fftw3.h>
#include <sstream>

#include "defines.hpp"
#include "PhaseSpace.hpp"
#include "Ruler.hpp"
#include "Impedance.hpp"
#include "CL/OpenCLHandler.hpp"

#include "IO/Display.hpp"

namespace vfps {

class ElectricField
{
public:
    /**
     * @brief ElectricField minimal constructor, will not offer wake function
     * @param phasespace this electric field is assigned to
     * @param impedance to use for electric field calculation
     * @param wakescaling scaling of wakepotential
     *        (As being part of fourier transform,
     *         delta t and delta f will be automatically taken into account.)
     *
     * Use other constructors when you want to use wake function or wake field.
     */
    ElectricField(PhaseSpace* ps,
                  const Impedance* impedance,
                  const double wakescalining=0.0);

    /**
     * @brief ElectricField constructor for use of wake potential
     * @param ps
     * @param impedance
     * @param Ib
     * @param E0
     * @param sigmaE
     * @param dt
     * @param rbend
     */
    ElectricField(PhaseSpace* ps,
                  const Impedance* impedance, const double Ib, const double E0,
                  const double sigmaE, const double dt, const double rbend);

    /**
     * @brief ElectricField (unmaintained) constructor for use of wake function
     * @param phasespace this electric field is assigned to
     * @param impedance to use for electric field calculation
     * @param Ib bunch current [A]
     * @param bl natural rms bunch length [m]
     * @param E0 beam energy [eV]
     * @param sigmaE normalized energy spread [1]
     * @param fs synchrotron frequency [Hz]
     * @param frev revolution frequency [Hz]
     * @param dt time step [s]
     * @param rbend bending radius [m]
     * @param nmax
     */
    ElectricField(PhaseSpace* ps, const Impedance* impedance,
                  const double Ib, const double E0,
                  const double sigmaE, const double fs,
                  const double dt, const double rbend,
                  const size_t nmax=0);

    ~ElectricField();

    inline csrpower_t* getData() const
        { return _csrspectrum; }

    inline const Impedance* getImpedance() const
        { return _impedance; }

    inline size_t getNMax() const
        { return _nmax; }

    inline const Ruler<meshaxis_t>* getFreqRuler() const
        { return &_axis_freq; }

    inline const Ruler<meshaxis_t>* getWakeRuler() const
        { return &_axis_wake; }

    csrpower_t* updateCSRSpectrum();

    meshaxis_t* getWakefunction() const
        { return _wakefunction; }

    /**
     * @brief wakePotential
     * @return
     *
     * @todo: Handling of negative frequencies in the formfactor
     * @todo: Correct scaling
     */
    meshaxis_t* wakePotential();

private:
    enum class fft_direction : uint_fast8_t {
        forward, backward
    };

    fftwf_plan prepareFFT(size_t n, csrpower_t* in, impedance_t* out);

    fftwf_plan prepareFFT(size_t n, impedance_t* in, impedance_t* out,
                          fft_direction direction);

private:
    const size_t _nmax;

    const size_t _bpmeshcells;

    const Ruler<meshaxis_t> _axis_freq;

    const Ruler<meshaxis_t> _axis_wake;

    PhaseSpace* _phasespace;

    csrpower_t* _csrspectrum;

    const Impedance* _impedance;

    integral_t* _bp_padded;

    integral_t* _bp_padded_fftw;

    #ifdef INOVESA_USE_CL
    cl::Buffer _bp_padded_clfft;
    #endif // INOVESA_USE_CL

    impedance_t* _formfactor;

    fftwf_complex* _formfactor_fftw;

    #ifdef INOVESA_USE_CL
    cl::Buffer _formfactor_clfft;
    #endif // INOVESA_USE_CL

    fftwf_plan _ffttw_bunchprofile;

    #ifdef INOVESA_USE_CL
    clfftPlanHandle _clfft_bunchprofile;
    #endif // INOVESA_USE_CL

    meshaxis_t* _wakefunction;

    impedance_t* _wakelosses;

    fftwf_complex* _wakelosses_fftw;

    #ifdef INOVESA_USE_CL
    cl::Buffer _wakelosses_clfft;
    #endif // INOVESA_USE_CL

    impedance_t* _wakepotential_complex;

    fftwf_complex* _wakepotential_fftw;

    #ifdef INOVESA_USE_CL
    cl::Buffer _wakepotential_clfft;
    #endif // INOVESA_USE_CL

    meshaxis_t* _wakepotential;

    fftwf_plan _fftwt_wakelosses;

    #ifdef INOVESA_USE_CL
    clfftPlanHandle _clfft_wakelosses;
    #endif // INOVESA_USE_CL

    const double _wakescaling;
};

} // namespace vfps

#endif // ELECTRICFIELD_HPP
