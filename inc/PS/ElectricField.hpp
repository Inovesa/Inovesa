/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2017: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2017: Karlsruhe Institute of Technology                 *
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
#include <memory>
#include <sstream>

#include "defines.hpp"
#include "PS/PhaseSpace.hpp"
#include "PS/Ruler.hpp"
#include "Z/Impedance.hpp"
#include "CL/OpenCLHandler.hpp"

#include "IO/Display.hpp"

namespace vfps {

typedef fftwf_plan fft_plan;
typedef fftwf_complex fft_complex;

class ElectricField
{
public:
    ElectricField()=delete;

    /**
     * @brief ElectricField minimal constructor, will not offer wake function
     * @param phasespace this electric field is assigned to
     * @param impedance to use for electric field calculation
     * @param revolutionpart
     * @param wakescaling scaling of wakepotential
     *        (As being part of fourier transform,
     *         delta t and delta f will be automatically taken into account.
     *         Also factor of 2 is applied to fix (Z*Rho)(k<0)=0)
     *
     * Use other constructors when you want to use wake function or potential.
     */
    ElectricField(std::shared_ptr<PhaseSpace> ps,
                  const std::shared_ptr<Impedance> impedance,
                  const double revolutionpart = 1,
                  const meshaxis_t wakescalining=0.0);

    /**
     * @brief ElectricField constructor for use of wake potential
     * @param phasespace this electric field is assigned to
     * @param impedance to use for electric field calculation
     * @param Ib bunch current [A]
     * @param E0 beam energy [eV]
     * @param sigmaE normalized energy spread [1]
     * @param dt time step [s]
     *
     * @todo: check whether impedance's frequencies match
     *
     * Internally all physical quantities are used to calculate
     * the wakescalining for ElectricField(PhaseSpace*,Impedance*,meshaxis_t).
     * We have factors of:
     *   Q_b = Ib/f0 (wake potential is calculated for normalized charge)
     *   dt*f0 (fraction of one revolution, impedance is for one revolution)
     *   c/ps->getScale(0) (charge/sigma_z -> current)
     *   1/(ps->getDelta(1)*sigmaE*E0) (eV -> pixels)
     */
    ElectricField(std::shared_ptr<PhaseSpace> ps,
                  std::shared_ptr<Impedance> impedance,
                  const double revolutionpart,
                  const double Ib, const double E0,
                  const double sigmaE, const double dt);

    /**
     * @brief ElectricField (unmaintained) constructor for use of wake function
     * @param phasespace this electric field is assigned to
     * @param impedance to use for electric field calculation
     * @param Ib bunch current [A]
     * @param bl natural rms bunch length [m]
     * @param E0 beam energy [eV]
     * @param sigmaE normalized energy spread [1]
     * @param frev revolution frequency [Hz]
     * @param dt time step [s]
     * @param fs synchrotron frequency [Hz]
     * @param nmax
     */
    ElectricField(std::shared_ptr<PhaseSpace> ps,
                  std::shared_ptr<Impedance> impedance,
                  const double Ib, const double E0,
                  const double sigmaE, const double dt, const double rbend,
                  const double fs, const size_t nmax);

    ~ElectricField();

    inline csrpower_t getCSRPower() const
        { return _csrintensity; }

    inline csrpower_t* getCSRSpectrum() const
        { return _csrspectrum; }

    inline const std::shared_ptr<Impedance> getImpedance() const
        { return _impedance; }

    inline size_t getNMax() const
        { return _nmax; }

    inline const Ruler<meshaxis_t>* getFreqRuler() const
        { return &_axis_freq; }

    inline const Ruler<meshaxis_t>* getWakeRuler() const
        { return &_axis_wake; }

    /**
     * @brief updateCSR
     * @param cutoff
     * @return CSR spectrum (getNMax() points)
     *
     * @todo: Use OpenCL for power calculation
     *
     * relies on an up-t date PhaseSpace::_projection[0]
     */
    csrpower_t* updateCSR(frequency_t cutoff);

    meshaxis_t* getWakefunction() const
        { return _wakefunction; }

    /**
     * @brief wakePotential
     * @return
     *
     * @todo: Handling of negative frequencies in the formfactor
     * @todo: Correct scaling
     *
     * relies on an up-t date PhaseSpace::_projection[axis]
     */
    meshaxis_t* wakePotential();

    #ifdef INOVESA_USE_CL
    void syncCLMem(clCopyDirection dir);
    #endif

public:
    const double volts;

private: // wrappers for FFTW
    enum class fft_direction : uint_fast8_t {
        forward, backward
    };

    fft_complex* fft_alloc_complex(size_t n);
    integral_t* fft_alloc_real(size_t n);

    /**
     * @brief fft_cleanup to be implemented
     */
    void fft_cleanup(){}

    inline void fft_destroy_plan(fftw_plan plan)
        { fftw_destroy_plan(plan); }
    inline void fft_destroy_plan(fftwf_plan plan)
        { fftwf_destroy_plan(plan); }

    inline void fft_execute(const fftw_plan plan)
        { fftw_execute(plan); }
    inline void fft_execute(const fftwf_plan plan)
        { fftwf_execute(plan); }

    inline void fft_free(double* addr)
        { fftw_free(addr); }
    inline void fft_free(float* addr)
        { fftwf_free(addr); }

    inline void fft_free(fftw_complex* addr)
        { fftw_free(addr); }
    inline void fft_free(fftwf_complex* addr)
        { fftwf_free(addr); }

    inline fftw_plan prepareFFT(size_t n, double* in, std::complex<double>* out)
        {return prepareFFT(n,in, reinterpret_cast<fftw_complex*>(out)); }
    fftw_plan prepareFFT(size_t n, double* in, fftw_complex* out);
    inline fftwf_plan prepareFFT(size_t n, float* in, std::complex<float>* out)
        {return prepareFFT(n,in, reinterpret_cast<fftwf_complex*>(out)); }
    fftwf_plan prepareFFT(size_t n, float* in, fftwf_complex* out);


    inline fftw_plan prepareFFT(size_t n, std::complex<double>* in,
                                std::complex<double>* out,
                                fft_direction direction)
        { return prepareFFT(n,reinterpret_cast<fftw_complex*>(in),
                            reinterpret_cast<fftw_complex*>(out),direction); }
    fftw_plan prepareFFT(size_t n, fftw_complex *in,
                         fftw_complex *out,
                         fft_direction direction);
    inline fftwf_plan prepareFFT(size_t n, std::complex<float>* in,
                                std::complex<float>* out,
                                fft_direction direction)
        { return prepareFFT(n,reinterpret_cast<fftwf_complex*>(in),
                            reinterpret_cast<fftwf_complex*>(out),direction); }
    fftwf_plan prepareFFT(size_t n, fftwf_complex* in,
                          fftwf_complex* out,
                          fft_direction direction);

private:
    const size_t _nmax;

    const uint32_t _bpmeshcells;

    const Ruler<meshaxis_t> _axis_freq;

    const Ruler<meshaxis_t> _axis_wake;

    std::shared_ptr<PhaseSpace> _phasespace;

    csrpower_t _csrintensity;

    csrpower_t* _csrspectrum;

    const std::shared_ptr<Impedance> _impedance;

    integral_t* _bp_padded;

    integral_t* _bp_padded_fft;

    #ifdef INOVESA_USE_CL
    cl::Buffer _bp_padded_buf;

    cl::Program _clProgPadBP;
    cl::Kernel _clKernPadBP;
    #endif // INOVESA_USE_CL

    impedance_t* _formfactor;

    fft_complex* _formfactor_fft;

    #ifdef INOVESA_USE_CL
    cl::Buffer _formfactor_buf;

    cl::Program _clProgWakelosses;
    cl::Kernel _clKernWakelosses;
    #endif // INOVESA_USE_CL

    fft_plan _fft_bunchprofile;

    #ifdef INOVESA_USE_CLFFT
    clfftPlanHandle _clfft_bunchprofile;
    #endif // INOVESA_USE_CLFFT

    meshaxis_t* _wakefunction;

    impedance_t* _wakelosses;

    fft_complex* _wakelosses_fft;

    #ifdef INOVESA_USE_CLFFT
    cl::Buffer _wakelosses_buf;
    #endif // INOVESA_USE_CLFFT

    impedance_t* _wakepotential_complex;

    fft_complex* _wakepotential_fft;

    #ifdef INOVESA_USE_CL
public:
    cl::Buffer _wakepotential_buf;

private:
    // non-interleaved internal data format might be usefull
    cl::Buffer _wakepotential_complex_buf;

    cl::Program _clProgScaleWP;
    cl::Kernel _clKernScaleWP;

    #endif // INOVESA_USE_CL

    meshaxis_t* _wakepotential;

    fft_plan _fft_wakelosses;

    #ifdef INOVESA_USE_CLFFT
    clfftPlanHandle _clfft_wakelosses;
    #endif // INOVESA_USE_CLFFT

    const meshdata_t _wakescaling;
};

} // namespace vfps

#endif // ELECTRICFIELD_HPP
