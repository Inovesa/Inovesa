// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include <algorithm>
#include <boost/multi_array.hpp>
#include <memory>
#include <sstream>

#include "defines.hpp"
#include "FFTWWrapper.hpp"
#include "PS/PhaseSpace.hpp"
#include "PS/Ruler.hpp"
#include "Z/Impedance.hpp"
#include "CL/OpenCLHandler.hpp"

#include "IO/Display.hpp"

namespace vfps {

class ElectricField
{
public:
    ElectricField()=delete;

    /**
     * @brief ElectricField minimal constructor, will not offer wake function
     * @param ps this electric field is assigned to
     * @param impedance to use for electric field calculation
     * @param bucketnumbers of buckets contained in ps
     * @param spacing_bins number of grid points between bunch centers
     *        (including grind points covered by phase space grid)
     * @param f_rev revolution frequency
     * @param revolutionpart part of one revolution covered by one time step
     * @param wakescaling scaling of wakepotential
     *        (As being part of fourier transform,
     *         delta t and delta f will be automatically taken into account.)
     *
     * Use other constructors when you want to use wake function or potential.
     */
    ElectricField(std::shared_ptr<PhaseSpace> ps,
                  const std::shared_ptr<Impedance> impedance,
                  const std::vector<uint32_t> &bucketnumbers,
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

    /**
     * @brief padBunchProfiles copies bunch profiles to have correct padding

     */
    void padBunchProfiles();

    inline integral_t* getPaddedBunchProfiles() const
        { return _bp_padded; }

    inline const auto& getWakePotentials() const
        { return _wakepotential; }

    inline meshaxis_t* getPaddedWakePotential() const
        { return _wakepotential_padded; }


    #if INOVESA_USE_OPENCL == 1
    void syncCLMem(OCLH::clCopyDirection dir);
    #endif // INOVESA_USE_OPENCL == 1

public:
    const double volts;

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
    const double factor4WattPerHertz;

    /**
     * @brief factor4Watts includes factor 2 to use positive frequencies only
     */
    const double factor4Watts;

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

    fft::complex* _formfactor_fft;

    oclhptr_t _oclh;

    #if INOVESA_USE_OPENCL == 1
    cl::Buffer _formfactor_buf;

    cl::Program _clProgWakelosses;
    cl::Kernel _clKernWakelosses;
    #endif // INOVESA_USE_OPENCL

    fft::plan _fft_bunchprofile;

    #if INOVESA_USE_CLFFT == 1
    clfftPlanHandle _clfft_bunchprofile;
    #endif // INOVESA_USE_CLFFT

    meshaxis_t* _wakefunction;

    impedance_t* _wakelosses;

    fft::complex* _wakelosses_fft;

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

    fft::plan _fft_wakelosses;

    #if INOVESA_USE_CLFFT == 1
    clfftPlanHandle _clfft_wakelosses;
    #endif // INOVESA_USE_CLFFT

    const meshdata_t _wakescaling;

public:
    /**
     * @brief getWakeScaling getter function for scaling factor
     */
    auto getWakeScaling() const
        { return _wakescaling; }
};

} // namespace vfps
