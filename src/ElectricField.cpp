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
                                   const size_t nmax) :
    _nmax(nmax > 0 ? nmax : impedance->maxN()),
    _axis(Ruler<meshaxis_t>(_nmax,0,meshaxis_t(1)/_nmax)),
    _phasespace(phasespace),
    _bpmeshcells(phasespace->nMeshCells(0)),
    _csrspectrum(new csrpower_t[_nmax]),
    _impedance(impedance),
    _spaceinfo(phasespace->getRuler(0)),
    _wakefunction(nullptr)
{
    _bp_padded_fftw = fftwf_alloc_real(2*_nmax);
    _bp_padded = reinterpret_cast<meshdata_t*>(_bp_padded_fftw);

    //zero-padding
    std::fill_n(&_bp_padded[_bpmeshcells],2*_nmax-_bpmeshcells,integral_t(0));

    _bp_fourier_fftw = fftwf_alloc_complex(2*_nmax);
    _bp_fourier = reinterpret_cast<impedance_t*>(_bp_fourier_fftw);

    _ft_bunchprofile = prepareFFT(2*_nmax,_bp_padded,_bp_fourier);
}

vfps::ElectricField::ElectricField(PhaseSpace* ps,
                                   const Impedance* impedance,
                                   const double Ib, const double bl,
                                   const double E0, const double sigmaE,
                                   const double fs, const double frev,
                                   const double dt, const double rbend,
                                   const size_t nmax) :
        ElectricField(ps,impedance,nmax)
{
    size_t wakenmax = std::round(physcons::c/(2*M_PI*frev*ps->getDelta(0)*bl));
    _wakefunction = new meshaxis_t[2*_bpmeshcells];
    fftwf_complex* z_fftw = fftwf_alloc_complex(wakenmax);
    fftwf_complex* zcsrf_fftw = fftwf_alloc_complex(wakenmax);
    fftwf_complex* zcsrb_fftw = fftwf_alloc_complex(wakenmax); //for wake
    impedance_t* z = reinterpret_cast<impedance_t*>(z_fftw);
    impedance_t* zcsrf = reinterpret_cast<impedance_t*>(zcsrf_fftw);
    impedance_t* zcsrb = reinterpret_cast<impedance_t*>(zcsrb_fftw);

    /* Marit's original code names eq1 in comment, but it uses eq2.
     * Patrik's calculations lead to eq3, but this does not work
     * right now, Patrik's result is adjusted the way Marit changed her result
     *
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
     *
     *  eq3:
     *  const double g    = -Ic * phaseSpace.getDelta<1>() / M_PI
     *                    * (deltat*omega0) * exp(sigma_z/R)
     *                    = -Ib * E0[eV] / (2*M_PI*f_s*sigma_delta)
     *                    * phaseSpace.getDelta<1>()
     *                    * (deltat*omega0 / M_PI) * exp(sigma_z/R)
     *                    = -Ib * E0[eV] / (2*M_PI*f_s*sigma_delta)
     *                    * phaseSpace.getDelta<1>()
     *                    * deltat*frev * exp(sigma_z/R)
     *                    = - Ib*E0*ps->getDelta(1)
     *                    * dt*frev*std::exp(bl/rbend)
     *                    / (2*M_PI*fs*sigmaE);
     */
     const double g = - Ib*physcons::c*ps->getDelta(1)*dt*std::exp(bl/rbend)
                    / (2*M_PI*fs*sigmaE*E0);


    std::copy_n(_impedance->data(),std::min(wakenmax,_impedance->maxN()),z);
    if (_impedance->maxN() < wakenmax) {
        std::stringstream wavenumbers;
        wavenumbers << "(Known: n=" <<_impedance->maxN()
                    << ", needed: N=" << wakenmax << ")";
        Display::printText("Warning: Unknown impedance for high wavenumbers. "
                           +wavenumbers.str());
        std::fill_n(&z[_impedance->maxN()],wakenmax-_impedance->maxN(),
                    impedance_t(0));
    }

    fftwf_plan p3 = prepareFFT( wakenmax, z, zcsrf, fft_direction::forward );
    fftwf_plan p4 = prepareFFT( wakenmax, z, zcsrb, fft_direction::backward);

    fftwf_execute(p3);
    fftwf_destroy_plan(p3);
    fftwf_execute(p4);
    fftwf_destroy_plan(p4);

    for (size_t i=0; i< _bpmeshcells; i++) {
        _wakefunction[i             ] = g * zcsrf[_bpmeshcells-i].real();
        _wakefunction[i+_bpmeshcells] = g * zcsrb[i             ].real();
    }
    fftwf_free(z_fftw);
    fftwf_free(zcsrf_fftw);
    fftwf_free(zcsrb_fftw);
}

vfps::ElectricField::~ElectricField()
{
    delete [] _csrspectrum;
    delete [] _wakefunction;

    fftwf_free(_bp_padded_fftw);
    fftwf_free(_bp_fourier_fftw);
    fftwf_destroy_plan(_ft_bunchprofile);
    fftwf_cleanup();
}

vfps::csrpower_t* vfps::ElectricField::updateCSRSpectrum()
{
    std::copy_n(_phasespace->projectionToX(),
                _spaceinfo->steps(),
                _bp_padded);

    //FFT charge density
    fftwf_execute(_ft_bunchprofile);

    for (unsigned int i=0; i<_nmax; i++) {
        // norm = squared magnitude
        _csrspectrum[i] = ((*_impedance)[i]).real()*std::norm(_bp_fourier[i]);
    }

    return _csrspectrum;
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
