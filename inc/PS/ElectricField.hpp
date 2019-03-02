// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include <algorithm>
#include <boost/multi_array.hpp>
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
     *
     * Use other constructors when you want to use wake function or potential.
     */
    ElectricField(std::shared_ptr<PhaseSpace> ps,
                  const std::shared_ptr<Impedance> impedance,
                  const std::vector<uint32_t> &bucketnumber,
                  const meshindex_t spacing_bins,
                  oclhptr_t oclh,
                  const double f_rev,
                  const meshaxis_t revolutionpart = 1,
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
    ElectricField( std::shared_ptr<PhaseSpace> ps
                 , std::shared_ptr<Impedance> impedance
                 , const std::vector<uint32_t> &bucketnumber
                 , const meshindex_t spacing_bins
                 , oclhptr_t oclh
                 , const double f_rev
                 , const double revolutionpart
                 , const double Ib, const double E0
                 , const double sigmaE, const double dt
                 );

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
    ElectricField(std::shared_ptr<PhaseSpace> ps
                 , std::shared_ptr<Impedance> impedance
                 , const std::vector<uint32_t> &bucketnumber
                 , const meshindex_t spacing_bins
                 , oclhptr_t oclh
                 , const double f_rev
                 , const double Ib, const double E0
                 , const double sigmaE, const double dt, const double rbend
                 , const double fs, const size_t nmax
                 );

    ~ElectricField() noexcept;

    inline const csrpower_t* getCSRPower() const
        { return _csrintensity.data(); }

    inline const csrpower_t* getCSRSpectrum() const
        { return _csrspectrum.data(); }

    inline const csrpower_t* getISRSpectrum() const
        { return _isrspectrum.data(); }

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
    const csrpower_t* updateCSR(const frequency_t cutoff);

    const std::vector<uint32_t> &getBuckets() const
        { return _bucket; }

    meshaxis_t* getWakefunction() const
        { return _wakefunction; }

    /**
     * @brief wakePotential updates wake potential
     * @return
     *
     * @todo: Handling of negative frequencies in the formfactor
     *
     * relies on an up to date PhaseSpace::_projection[axis]
     */
    meshaxis_t* wakePotential();

    inline integral_t* getPaddedProfile() const
        { return _bp_padded; }

    inline meshaxis_t* getPaddedWakepotential() const
        { return _wakepotential_padded; }


    #if INOVESA_USE_OPENCL == 1
    void syncCLMem(OCLH::clCopyDirection dir);
    #endif // INOVESA_USE_OPENCL == 1

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

    /**
     * @brief prepareFFT real to complex (Hermitian) FFT, "forward"
     * @param n
     * @param in
     * @param out
     * @return
     */
    fftwf_plan prepareFFT(size_t n, float* in, fftwf_complex* out);

    inline fftwf_plan prepareFFT(size_t n, std::complex<float>* in, float* out)
        {return prepareFFT(n,reinterpret_cast<fftwf_complex*>(in), out); }

    /**
     * @brief prepareFFT Hermitian (complex) to real FFT, "backward"
     * @param n
     * @param in
     * @param out
     * @return
     */
    fftwf_plan prepareFFT(size_t n, fftwf_complex* in, float* out);

    inline fftw_plan prepareFFT(size_t n, std::complex<double>* in,
                                std::complex<double>* out,
                                fft_direction direction)
        { return prepareFFT(n,reinterpret_cast<fftw_complex*>(in),
                            reinterpret_cast<fftw_complex*>(out),direction); }

    /**
     * @brief prepareFFT (unmaintained for C2C FFT)
     */
    fftw_plan prepareFFT(size_t n, fftw_complex *in,
                         fftw_complex *out,
                         fft_direction direction);

    inline fftwf_plan prepareFFT(size_t n, std::complex<float>* in,
                                std::complex<float>* out,
                                fft_direction direction)
        { return prepareFFT(n,reinterpret_cast<fftwf_complex*>(in),
                            reinterpret_cast<fftwf_complex*>(out),direction); }

    /**
     * @brief prepareFFT (unmaintained for C2C FFT)
     */
    fftwf_plan prepareFFT(size_t n, fftwf_complex* in,
                          fftwf_complex* out,
                          fft_direction direction);

private:
    const uint32_t _nbunches;

    const std::vector<uint32_t> _bucket;

    const size_t _nmax;

    const size_t _spacing_bins;

    const Ruler<meshaxis_t> _axis_freq;

    const Ruler<meshaxis_t> _axis_wake;

    std::shared_ptr<PhaseSpace> _phasespace;

    /// factor to normalize form factor
    const csrpower_t _formfactorrenorm;

public:
    /**
     * @brief OhmsPerHertz includes factor 2 to use positive frequencies only
     */
    const csrpower_t factor4WattPerHertz;

    /**
     * @brief factor4Watts includes factor 2 to use positive frequencies only
     */
    const csrpower_t factor4Watts;

private:


    /**
     * @brief _csrintensity dimensions: bunch
     */
    boost::multi_array<csrpower_t,1> _csrintensity;

    /**
     * @brief _csrspectrum dimensions: bunch, frequency
     */
    boost::multi_array<csrpower_t,2> _csrspectrum;

    /**
     * @brief _isrspectrum dimensions: bunch, frequency
     */
    boost::multi_array<csrpower_t,2> _isrspectrum;

    const std::shared_ptr<Impedance> _impedance;

    integral_t* _bp_padded;

    integral_t* _bp_padded_fft;

    #if INOVESA_USE_OPENCL == 1
    cl::Buffer _bp_padded_buf;
    #endif // INOVESA_USE_OPENCL

    impedance_t* _formfactor;

    fft_complex* _formfactor_fft;

    oclhptr_t _oclh;

    #if INOVESA_USE_OPENCL == 1
    cl::Buffer _formfactor_buf;

    cl::Program _clProgWakelosses;
    cl::Kernel _clKernWakelosses;
    #endif // INOVESA_USE_OPENCL

    fft_plan _fft_bunchprofile;

    #if INOVESA_USE_CLFFT == 1
    clfftPlanHandle _clfft_bunchprofile;
    #endif // INOVESA_USE_CLFFT

    meshaxis_t* _wakefunction;

    impedance_t* _wakelosses;

    fft_complex* _wakelosses_fft;

    #if INOVESA_USE_CLFFT == 1
    cl::Buffer _wakelosses_buf;
    #endif // INOVESA_USE_CLFFT

    /**
     * @brief _wakepotential_complex wake potential of size _nmax
     *
     * Actually, this is a real value.
     * Implement usage of C2R FFT to use that fact.
     */
    meshaxis_t* _wakepotential_padded;

    #if INOVESA_USE_OPENCL == 1
public:
    #if INOVESA_USE_OPENGL == 1
    cl_GLuint wakepotential_glbuf;
    #endif // INOVESA_USE_OPENGL

    cl::Buffer wakepotential_clbuf;

private:
    // @todo: non-interleaved internal data format might be usefull
    cl::Buffer _wakepotential_padded_buf;

    cl::Program _clProgScaleWP;
    cl::Kernel _clKernScaleWP;

    #endif // INOVESA_USE_OPENCL

    boost::multi_array<meshaxis_t,2> _wakepotential;

    fft_plan _fft_wakelosses;

    #if INOVESA_USE_CLFFT == 1
    clfftPlanHandle _clfft_wakelosses;
    #endif // INOVESA_USE_CLFFT

    const meshdata_t _wakescaling;
};

} // namespace vfps
