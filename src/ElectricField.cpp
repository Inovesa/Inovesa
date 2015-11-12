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

vfps::ElectricField::ElectricField(PhaseSpace* phasespace,
                                   const Impedance* impedance,
                                   const bool wakepot) :
    _nmax(impedance->maxN()),
    _bpmeshcells(phasespace->nMeshCells(0)),
    _axis_freq(Ruler<meshaxis_t>(_nmax,0,meshaxis_t(1)/_nmax)),
    // _axis_wake[_bpmeshcells] will be 0
    _axis_wake(Ruler<meshaxis_t>(2*_bpmeshcells,
                                 -phasespace->getDelta(0)*_bpmeshcells,
                                  phasespace->getDelta(0)*(_bpmeshcells-1),
                                 phasespace->getScale(0))),
    _phasespace(phasespace),
    _csrspectrum(new csrpower_t[_nmax]),
    _impedance(impedance),
    _wakefunction(nullptr),
    _wakepotential(wakepot?new meshaxis_t[_bpmeshcells]:nullptr)
{
    _bp_padded_fftw = fftwf_alloc_real(_nmax);
    _bp_padded = reinterpret_cast<meshdata_t*>(_bp_padded_fftw);
    std::fill_n(_bp_padded,_nmax,integral_t(0));

    _bp_fourier_fftw = fftwf_alloc_complex(_nmax);
    _bp_fourier = reinterpret_cast<impedance_t*>(_bp_fourier_fftw);
    std::fill_n(_bp_fourier,_nmax,integral_t(0));

    _ft_bunchprofile = prepareFFT(_nmax,_bp_padded,_bp_fourier);


    if (wakepot) {
        _wakelosses_fftw = fftwf_alloc_complex(_nmax);
        _wakepotential_fftw = fftwf_alloc_complex(_nmax);
    } else {
        _wakelosses_fftw = nullptr;
        _wakepotential_fftw = nullptr;
    }

    _wakelosses=reinterpret_cast<impedance_t*>(_wakelosses_fftw);
    _wakepotential_complex=reinterpret_cast<impedance_t*>(_wakepotential_fftw);

    if (wakepot) {
        _ft_wakelosses = prepareFFT(_nmax,_wakelosses,
                                    _wakepotential_complex,
                                    fft_direction::forward);
    } else {
        _ft_wakelosses = nullptr;
    }
}

vfps::ElectricField::ElectricField(PhaseSpace* ps, const Impedance* impedance,
                                   const double Ib, const double E0,
                                   const double sigmaE, const double fs,
                                   const double dt, const double rbend,
                                   const size_t nmax) :
        ElectricField(ps,impedance,nmax)
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

    fftwf_free(_bp_padded_fftw);
    fftwf_free(_bp_fourier_fftw);
    if(_wakelosses_fftw != nullptr) {
        fftwf_free(_wakelosses_fftw);
    }
    if(_wakepotential_fftw != nullptr) {
        fftwf_free(_wakepotential_fftw);
    }
    fftwf_destroy_plan(_ft_bunchprofile);
    if (_ft_wakelosses != nullptr) {
        fftwf_destroy_plan(_ft_wakelosses);
    }
    fftwf_cleanup();
}

vfps::csrpower_t* vfps::ElectricField::updateCSRSpectrum()
{
    // copy bunch profile so that negative times are at maximum bins
    vfps::projection_t* bp=_phasespace->projectionToX();
    std::copy_n(bp,_bpmeshcells/2,_bp_padded+_nmax-_bpmeshcells/2);
    std::copy_n(bp+_bpmeshcells/2,_bpmeshcells/2,_bp_padded);

    //FFT charge density
    fftwf_execute(_ft_bunchprofile);

    for (unsigned int i=0; i<_nmax; i++) {
        // norm = squared magnitude
        _csrspectrum[i] = ((*_impedance)[i]).real()*std::norm(_bp_fourier[i]);
    }

    return _csrspectrum;
}

vfps::meshaxis_t *vfps::ElectricField::wakePotential()
{
    // copy bunch profile so that negative times are at maximum bins
    vfps::projection_t* bp=_phasespace->projectionToX();
    std::copy_n(bp,_bpmeshcells/2,_bp_padded+_nmax-_bpmeshcells/2);
    std::copy_n(bp+_bpmeshcells/2,_bpmeshcells/2,_bp_padded);

    //FFT charge density
    fftwf_execute(_ft_bunchprofile);

    std::fill_n(_wakelosses,_nmax,0);
    for (unsigned int i=0; i<_nmax; i++) {
        _wakelosses[i]=((*_impedance)[i]*_bp_fourier[i]);
    }

    fftwf_execute(_ft_wakelosses);

    for (unsigned int i=0; i<_bpmeshcells/2; i++) {
        _wakepotential[_bpmeshcells/2+i]=_wakepotential_complex[i].real();
        _wakepotential[_bpmeshcells/2-i-1]
            =_wakepotential_complex[_nmax-i-1].real();
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
