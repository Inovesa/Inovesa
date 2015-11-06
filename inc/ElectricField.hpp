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

#include "IO/Display.hpp"

namespace vfps {

class ElectricField
{
public:
    /**
     * @brief ElectricField minimal constructor, will not offer wake function
     * @param phasespace this electric field is assigned to
     * @param impedance to use for electric field calculation
     */
    ElectricField(PhaseSpace* phasespace,
                  const Impedance* impedance,
                  const size_t nmax=0, const size_t padding=2,
                  const bool wakepot=false);

    /**
     * @brief ElectricField
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

    const size_t _padding;

    const size_t _bpmeshcells;

    const Ruler<meshaxis_t> _axis_freq;

    const Ruler<meshaxis_t> _axis_wake;

    PhaseSpace* _phasespace;

    csrpower_t* _csrspectrum;

    const Impedance* _impedance;

    integral_t* _bp_padded;

    integral_t* _bp_padded_fftw;

    impedance_t* _bp_fourier;

    fftwf_complex* _bp_fourier_fftw;

    fftwf_plan _ft_bunchprofile;

    meshaxis_t* _wakefunction;

    impedance_t* _wakelosses;

    fftwf_complex* _wakelosses_fftw;

    impedance_t* _wakepotential_complex;

    fftwf_complex* _wakepotential_fftw;

    meshaxis_t* _wakepotential;

    fftwf_plan _ft_wakelosses;
};

} // namespace vfps

#endif // ELECTRICFIELD_HPP
