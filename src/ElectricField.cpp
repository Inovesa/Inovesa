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

#include "ElectricField.hpp"

vfps::ElectricField::ElectricField(PhaseSpace* ps,
                                   const Impedance* impedance,
                                   double wakescalining) :
    _nmax(impedance->maxN()),
    _bpmeshcells(ps->nMeshCells(0)),
    _axis_freq(Ruler<frequency_t>(_nmax,0,1/(ps->getDelta(0)),0)),
    // _axis_wake[_bpmeshcells] will be at position 0
    _axis_wake(Ruler<meshaxis_t>(2*_bpmeshcells,
                                 -ps->getDelta(0)*_bpmeshcells,
                                  ps->getDelta(0)*(_bpmeshcells-1),
                                 ps->getScale(0))),
    _phasespace(ps),
    _csrspectrum(new csrpower_t[_nmax]),
    _impedance(impedance),
    _wakefunction(nullptr),
    _wakelosses(nullptr),
    _wakelosses_fftw(nullptr),
    _wakepotential_complex(nullptr),
    _wakepotential_fftw(nullptr),
    _wakepotential(wakescalining!=0.0?new meshaxis_t[_bpmeshcells]:nullptr),
    _fftwt_wakelosses(nullptr),
    _wakescaling(wakescalining*_axis_wake.delta()*_axis_freq.delta())
{
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        _bp_padded = new meshdata_t[_nmax];
        std::fill_n(_bp_padded,_nmax,0);
        _bp_padded_clfft = cl::Buffer(OCLH::context,
                                      CL_MEM_READ_WRITE,
                                      sizeof(impedance_t)*_nmax,_bp_padded);
        _formfactor = new impedance_t[_nmax];
        _formfactor_clfft = cl::Buffer(OCLH::context,
                                       CL_MEM_READ_WRITE,
                                       sizeof(impedance_t)*_nmax,_formfactor);
        clfftCreateDefaultPlan(&_clfft_bunchprofile,
                               OCLH::context(),CLFFT_1D,&_nmax);
        clfftSetPlanPrecision(_clfft_bunchprofile,CLFFT_SINGLE);
        clfftSetLayout(_clfft_bunchprofile, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED);
        clfftSetResultLocation(_clfft_bunchprofile, CLFFT_OUTOFPLACE);
        clfftBakePlan(_clfft_bunchprofile,1,&OCLH::queue(), nullptr, nullptr);

        // to be implemented
        // std::copy_n(bp,_bpmeshcells/2,_bp_padded+_nmax-_bpmeshcells/2);
        // std::copy_n(bp+_bpmeshcells/2,_bpmeshcells/2,_bp_padded);
        std::string cl_code_padbp = R"(
            __kernel void pad_bp(__global float* bp_padded,
                                 const uint paddedsize,
                                 const uint bpmeshcells,
                                 const __global float* bp)
            {
                const uint g = get_global_id(0);
                const uint b = (g+bpmeshcells/2)%bpmeshcells;
                const uint p;
                bp_padded[p] = bp[b];
            }
            )";

        _clProgPadBP = OCLH::prepareCLProg(cl_code_padbp);
        _clKernPadBP = cl::Kernel(_clProgPadBP, "pad_bp");
        _clKernPadBP.setArg(0, _bp_padded_clfft);
        _clKernPadBP.setArg(1, _nmax);
        _clKernPadBP.setArg(2, _bpmeshcells);
        _clKernPadBP.setArg(3, _phasespace->projectionX_buf);
    } else
    #endif // INOVESA_USE_CL
    {
        _bp_padded_fftw = fftwf_alloc_real(_nmax);
        _bp_padded = reinterpret_cast<meshdata_t*>(_bp_padded_fftw);
        std::fill_n(_bp_padded,_nmax,integral_t(0));
        _formfactor_fftw = fftwf_alloc_complex(_nmax);
        _formfactor = reinterpret_cast<impedance_t*>(_formfactor_fftw);
        std::fill_n(_formfactor,_nmax,integral_t(0));
        _ffttw_bunchprofile = prepareFFT(_nmax,_bp_padded,_formfactor);
    }
}

vfps::ElectricField::ElectricField(vfps::PhaseSpace *ps,
                                   const vfps::Impedance *impedance,
                                   const double Ib, const double E0,
                                   const double sigmaE, const double dt,
                                   const double rbend) :
    ElectricField(ps,impedance,4*M_PI*rbend*Ib/physcons::c*dt/physcons::e
                               /(ps->getDelta(1)*sigmaE*E0)
                 )
{
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        _wakelosses = new impedance_t[_nmax];
        _wakelosses_clfft = cl::Buffer(OCLH::context,
                                      CL_MEM_READ_WRITE,
                                      sizeof(impedance_t)*_nmax,_wakelosses);
        _wakepotential = new meshaxis_t[_nmax];
        _wakepotential_clfft = cl::Buffer(OCLH::context,
                                          CL_MEM_READ_WRITE,
                                          sizeof(meshaxis_t)*_nmax,
                                          _wakepotential);
        clfftCreateDefaultPlan(&_clfft_wakelosses,
                               OCLH::context(),CLFFT_1D,&_nmax);
        clfftSetPlanPrecision(_clfft_wakelosses,CLFFT_SINGLE);
        clfftSetLayout(_clfft_wakelosses,CLFFT_HERMITIAN_INTERLEAVED,CLFFT_REAL);
        clfftSetResultLocation(_clfft_wakelosses, CLFFT_OUTOFPLACE);
        clfftBakePlan(_clfft_wakelosses,1,&OCLH::queue(), nullptr, nullptr);

        std::string cl_code_wakelosses = R"(
            __kernel void wakeloss(__global float* wakelosses,
                                   const float scaling,
                                   const __global float* impedance,
                                   const __global float* formfactor)
            {
                const uint n = get_global_id(0);
                wakelosses[n] = scaling*impedance[n]*formfactor[n];
            }
            )";

        _clProgPadBP = OCLH::prepareCLProg(cl_code_wakelosses);
        _clKernWakelosses = cl::Kernel(_clProgPadBP, "wakeloss");
        _clKernWakelosses.setArg(0, _wakelosses_clfft);
        _clKernWakelosses.setArg(1, _wakescaling);
        _clKernWakelosses.setArg(2, _impedance->data_buf);
        _clKernWakelosses.setArg(3, _formfactor_clfft);
    } else
    #endif // INOVESA_USE_CL
    {
        _wakelosses_fftw = fftwf_alloc_complex(_nmax);
        _wakepotential_fftw = fftwf_alloc_complex(_nmax);

        _wakelosses=reinterpret_cast<impedance_t*>(_wakelosses_fftw);
        _wakepotential_complex=reinterpret_cast<impedance_t*>(_wakepotential_fftw);
        _fftwt_wakelosses = prepareFFT(_nmax,_wakelosses,
                                    _wakepotential_complex,
                                  fft_direction::backward);
    }
}

// (unmaintained) constructor for use of wake function
vfps::ElectricField::ElectricField(PhaseSpace* ps, const Impedance* impedance,
                                   const double Ib, const double E0,
                                   const double sigmaE, const double fs,
                                   const double dt, const double rbend,
                                   const size_t nmax) :
        ElectricField(ps,impedance)
{
    _wakefunction = new meshaxis_t[2*_bpmeshcells];
    fftwf_complex* z_fftw = fftwf_alloc_complex(nmax);
    fftwf_complex* zcsrf_fftw = fftwf_alloc_complex(nmax);
    fftwf_complex* zcsrb_fftw = fftwf_alloc_complex(nmax); //for wake
    impedance_t* z = reinterpret_cast<impedance_t*>(z_fftw);
    impedance_t* zcsrf = reinterpret_cast<impedance_t*>(zcsrf_fftw);
    impedance_t* zcsrb = reinterpret_cast<impedance_t*>(zcsrb_fftw);

    /* Marit's original code names eq1 in comment, but it uses eq2.
     *
     *    eq1:
     *     const double g = -Ic * phaseSpace.getDelta<0>() / M_PI
     *                    * (deltat*omega0);
     *
     *    eq2:
     *    const double g  = -Ic * phaseSpace.getDelta<0>() / M_PI
     *                    * deltat * physcons::c * beta0 / R
     *
     * Marit's comment:
     * !!! omega0 is here a function of R !!!, deltat in Einheiten von 2*pi?
     */
     const double g = - Ib*physcons::c*ps->getDelta(1)*dt
                    / (2*M_PI*fs*sigmaE*E0)/(M_PI*rbend);


    std::copy_n(_impedance->data(),std::min(nmax,_impedance->maxN()),z);
    if (_impedance->maxN() < nmax) {
        std::stringstream wavenumbers;
        wavenumbers << "(Known: n=" <<_impedance->maxN()
                    << ", needed: N=" << nmax << ")";
        Display::printText("Warning: Unknown impedance for high wavenumbers. "
                           +wavenumbers.str());
        std::fill_n(&z[_impedance->maxN()],nmax-_impedance->maxN(),
                    impedance_t(0));
    }

    fftwf_plan p3 = prepareFFT( nmax, z, zcsrf, fft_direction::forward );
    fftwf_plan p4 = prepareFFT( nmax, z, zcsrb, fft_direction::backward);

    fftwf_execute(p3);
    fftwf_destroy_plan(p3);
    fftwf_execute(p4);
    fftwf_destroy_plan(p4);

    /* This method works like a DFT of Z with Z(-n) = Z*(n).
     *
     * the element _wakefunction[_bpmeshcells] represents the self interaction
     * set this element (q==0) to zero to make the function anti-semetric
     */
    _wakefunction[0] = 0;
    for (size_t i=0; i< _bpmeshcells; i++) {
        // zcsrf[0].real() == zcsrb[0].real(), see comment above
        _wakefunction[_bpmeshcells-i] = g * zcsrf[i].real();
        _wakefunction[_bpmeshcells+i] = g * zcsrb[i].real();
    }
    fftwf_free(z_fftw);
    fftwf_free(zcsrf_fftw);
    fftwf_free(zcsrb_fftw);
}

vfps::ElectricField::~ElectricField()
{
    delete [] _csrspectrum;
    delete [] _wakefunction;
    delete [] _wakepotential;

    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        delete [] _bp_padded;
        delete [] _formfactor;
        delete [] _wakepotential;
        clfftDestroyPlan(&_clfft_bunchprofile);
        clfftDestroyPlan(&_clfft_wakelosses);
    } else
    #endif // INOVESA_USE_CL
    {
        fftwf_free(_bp_padded_fftw);
        fftwf_free(_formfactor_fftw);
        if(_wakelosses_fftw != nullptr) {
            fftwf_free(_wakelosses_fftw);
        }
        if(_wakepotential_fftw != nullptr) {
            fftwf_free(_wakepotential_fftw);
        }
        fftwf_destroy_plan(_ffttw_bunchprofile);
        if (_fftwt_wakelosses != nullptr) {
            fftwf_destroy_plan(_fftwt_wakelosses);
        }
        fftwf_cleanup();
    }
}

vfps::csrpower_t* vfps::ElectricField::updateCSRSpectrum()
{
    _phasespace->updateXProjection();
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        clfftEnqueueTransform(_clfft_bunchprofile,CLFFT_FORWARD,1,&OCLH::queue(),
                          0,nullptr,nullptr,
                          &_bp_padded_clfft(),&_formfactor_clfft(),nullptr);
        OCLH::queue.enqueueBarrierWithWaitList();
    } else
    #endif // INOVESA_USE_CL
    {
        // copy bunch profile so that negative times are at maximum bins
        const vfps::projection_t* bp= _phasespace->getProjection(0);
        std::copy_n(bp,_bpmeshcells/2,_bp_padded+_nmax-_bpmeshcells/2);
        std::copy_n(bp+_bpmeshcells/2,_bpmeshcells/2,_bp_padded);
        //FFT charge density
        fftwf_execute(_ffttw_bunchprofile);

        for (unsigned int i=0; i<_nmax; i++) {
            // norm = squared magnitude
            _csrspectrum[i] = ((*_impedance)[i]).real()*std::norm(_formfactor[i]);
        }
    }

    return _csrspectrum;
}

vfps::meshaxis_t *vfps::ElectricField::wakePotential()
{
    _phasespace->updateXProjection();
    #ifdef INOVESA_USE_CL
    if (OCLH::active){
        OCLH::queue.enqueueNDRangeKernel( _clKernPadBP,cl::NullRange,
                                          cl::NDRange(_bpmeshcells));
        OCLH::queue.enqueueBarrierWithWaitList();
        clfftEnqueueTransform(_clfft_bunchprofile,CLFFT_FORWARD,1,&OCLH::queue(),
                          0,nullptr,nullptr,
                          &_bp_padded_clfft(),&_formfactor_clfft(),nullptr);
        OCLH::queue.enqueueBarrierWithWaitList();
        OCLH::queue.enqueueNDRangeKernel( _clKernWakelosses,cl::NullRange,
                                          cl::NDRange(_nmax));
        OCLH::queue.enqueueBarrierWithWaitList();
        clfftEnqueueTransform(_clfft_wakelosses,CLFFT_BACKWARD,1,&OCLH::queue(),
                          0,nullptr,nullptr,
                          &_wakelosses_clfft(),&_wakepotential_clfft(),nullptr);
        OCLH::queue.enqueueBarrierWithWaitList();
    } else
    #endif // INOVESA_USE_CL
    {
        // copy bunch profile so that negative times are at maximum bins
        const vfps::projection_t* bp= _phasespace->getProjection(0);
        std::copy_n(bp,_bpmeshcells/2,_bp_padded+_nmax-_bpmeshcells/2);
        std::copy_n(bp+_bpmeshcells/2,_bpmeshcells/2,_bp_padded);

        // Fourier transorm charge density
        // FFTW R2C only computes elements 0...n/2, and
        // sets second half of output array to 0.
        // This is because Y[n-i] = Y[i].
        // We will use this, and choose the wake losses
        // for negetive frequencies to be 0, equivalent to Z(-|f|)=0.
        fftwf_execute(_ffttw_bunchprofile);

        for (unsigned int i=0; i<_nmax/2; i++) {
            _wakelosses[i]= (*_impedance)[i] *_formfactor[i];
        }
        std::fill_n(_wakelosses+_nmax/2,_nmax/2,0);

        //Fourier transorm wakelosses
        fftwf_execute(_fftwt_wakelosses);

        for (size_t i=0; i<_bpmeshcells/2; i++) {
            _wakepotential[_bpmeshcells/2+i]
                = _wakepotential_complex[        i].real()*_wakescaling;
            _wakepotential[_bpmeshcells/2-1-i]
                = _wakepotential_complex[_nmax-1-i].real()*_wakescaling;
        }
    }

    return _wakepotential;
}

fftwf_plan vfps::ElectricField::prepareFFT( size_t n, csrpower_t* in,
                                            impedance_t* out)
{
    fftwf_plan plan = nullptr;

    std::stringstream wisdomfile;
    // get ready to save BunchCharge
    if (std::is_same<vfps::csrpower_t,float>::value) {
        wisdomfile << "wisdom_r2c32_" << n << ".fftw";
    } else if (std::is_same<vfps::csrpower_t,double>::value) {
        wisdomfile << "wisdom_r2c64_" << n << ".fftw";
    }
    // use wisdomfile, if it exists
    if (fftwf_import_wisdom_from_filename(wisdomfile.str().c_str()) != 0) {
        plan = fftwf_plan_dft_r2c_1d(n,in,reinterpret_cast<fftwf_complex*>(out),
                                    FFTW_WISDOM_ONLY|FFTW_PATIENT);
    }
    // if there was no wisdom (no or bad file), create some
    if (plan == nullptr) {
        plan = fftwf_plan_dft_r2c_1d(n,in,reinterpret_cast<fftwf_complex*>(out),
                                    FFTW_PATIENT);
        fftwf_export_wisdom_to_filename(wisdomfile.str().c_str());
        Display::printText("Created some wisdom at "+wisdomfile.str());
    }
    return plan;
}

fftwf_plan vfps::ElectricField::prepareFFT(size_t n, vfps::impedance_t* in,
                                    vfps::impedance_t* out,
                                    fft_direction direction)
{
    fftwf_plan plan = nullptr;

    char dir;
    int_fast8_t sign;
    if (direction == fft_direction::backward) {
        dir = 'b';
        sign = +1;
    } else {
        dir = 'f';
        sign = -1;
    }

    std::stringstream wisdomfile;
    // get ready to save BunchCharge
    if (std::is_same<vfps::csrpower_t,float>::value) {
        wisdomfile << "wisdom_c" << dir << "c32_" << n << ".fftw";
    } else {
        wisdomfile << "wisdom_c" << dir << "c64_" << n << ".fftw";
    }
    // use wisdomfile, if it exists
    if (fftwf_import_wisdom_from_filename(wisdomfile.str().c_str()) != 0) {
        plan = fftwf_plan_dft_1d(    n,
                                        reinterpret_cast<fftwf_complex*>(in),
                                        reinterpret_cast<fftwf_complex*>(out),
                                        sign,
                                        FFTW_WISDOM_ONLY|FFTW_PATIENT);
    }
    // if there was no wisdom (no or bad file), create some
    if (plan == nullptr) {
        plan = fftwf_plan_dft_1d(    n,
                                        reinterpret_cast<fftwf_complex*>(in),
                                        reinterpret_cast<fftwf_complex*>(out),
                                        sign,
                                        FFTW_PATIENT);
        fftwf_export_wisdom_to_filename(wisdomfile.str().c_str());
        Display::printText("Created some wisdom at "+wisdomfile.str());
    }
    return plan;

}
