// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#include "FFTWWrapper.hpp"

#include "IO/Display.hpp"
#include "IO/FSPath.hpp"

fft::complex* fft::fft_alloc_complex(size_t n)
{
    fft::complex* rv;
    if (std::is_same<vfps::csrpower_t,float>::value) {
        rv = reinterpret_cast<complex*>(fftwf_alloc_complex(n));
    } else if (std::is_same<vfps::csrpower_t,double>::value) {
        rv = reinterpret_cast<complex*>(fftw_alloc_complex(n));
    }
    // initialize as 2*n long real array
    std::fill_n(reinterpret_cast<vfps::csrpower_t*>(rv),2*n,0);
    return rv;
}

vfps::integral_t* fft::fft_alloc_real(size_t n)
{
    vfps::integral_t* rv;
    if (std::is_same<vfps::csrpower_t,float>::value) {
        rv = reinterpret_cast<vfps::integral_t*>(fftwf_alloc_real(n));
    } else if (std::is_same<vfps::csrpower_t,double>::value) {
        rv = reinterpret_cast<vfps::integral_t*>(fftw_alloc_real(n));
    }
    std::fill_n(rv,n,vfps::integral_t(0));
    return rv;
}

fftw_plan fft::prepareFFT(size_t n, double* in,
                                          fftw_complex* out)
{
    fftw_plan plan = nullptr;

    std::stringstream wisdomfname;
    wisdomfname << "wisdom_r2c64_" << n << ".fftw";

    vfps::FSPath wisdompath(vfps::FSPath::datapath());
    wisdompath.append("fftwisdom/"+wisdomfname.str());

    // use wisdomfile, if it exists
    if (fftw_import_wisdom_from_filename(wisdompath.c_str()) != 0) {
        plan = fftw_plan_dft_r2c_1d(n,in,out,FFTW_WISDOM_ONLY|FFTW_PATIENT);
    }
    // if there was no wisdom (no or bad file), create some
    if (plan == nullptr) {
        plan = fftw_plan_dft_r2c_1d(n,in,out,FFTW_PATIENT);
        fftw_export_wisdom_to_filename(wisdompath.c_str());
        vfps::Display::printText("Created some wisdom at "+wisdompath.str());
    }
    return plan;
}


fftwf_plan fft::prepareFFT( size_t n, float *in,
                            fftwf_complex*out)
{
    fftwf_plan plan = nullptr;

    std::stringstream wisdomfname;
    wisdomfname << "wisdom_r2c32_" << n << ".fftw";

    vfps::FSPath wisdompath(vfps::FSPath::datapath());
    wisdompath.append("fftwisdom/"+wisdomfname.str());

    // use wisdomfile, if it exists
    if (fftwf_import_wisdom_from_filename(wisdompath.c_str()) != 0) {
        plan = fftwf_plan_dft_r2c_1d(n,in,out,FFTW_WISDOM_ONLY|FFTW_PATIENT);
    }
    // if there was no wisdom (no or bad file), create some
    if (plan == nullptr) {
        plan = fftwf_plan_dft_r2c_1d(n,in,out,FFTW_PATIENT);
        fftwf_export_wisdom_to_filename(wisdompath.c_str());
        vfps::Display::printText("Created some wisdom at "+wisdompath.str());
    }
    return plan;
}

fftwf_plan fft::prepareFFT(size_t n, fftwf_complex *in,
                                           float *out)
{
    fftwf_plan plan = nullptr;

    std::stringstream wisdomfname;
    wisdomfname << "wisdom_c2r32_" << n << ".fftw";

    vfps::FSPath wisdompath(vfps::FSPath::datapath());
    wisdompath.append("fftwisdom/"+wisdomfname.str());

    // use wisdomfile, if it exists
    if (fftwf_import_wisdom_from_filename(wisdompath.c_str()) != 0) {
        plan = fftwf_plan_dft_c2r_1d(n,in,out,FFTW_WISDOM_ONLY|FFTW_PATIENT);
    }
    // if there was no wisdom (no or bad file), create some
    if (plan == nullptr) {
        plan = fftwf_plan_dft_c2r_1d(n,in,out,FFTW_PATIENT);
        fftwf_export_wisdom_to_filename(wisdompath.c_str());
        vfps::Display::printText("Created some wisdom at "+wisdompath.str());
    }
    return plan;
}

fftw_plan fft::prepareFFT( size_t n,
                           fftw_complex* in,
                           fftw_complex* out,
                           fft_direction direction)
{
    fftw_plan plan = nullptr;

    char dir;
    int_fast8_t sign;
    if (direction == fft_direction::backward) {
        dir = 'b';
        sign = +1;
    } else {
        dir = 'f';
        sign = -1;
    }

    std::stringstream wisdomfname;
    // find filename for wisdom
    wisdomfname << "wisdom_c" << dir << "c64_" << n << ".fftw";

    vfps::FSPath wisdompath(vfps::FSPath::datapath());
    wisdompath.append("fftwisdom/"+wisdomfname.str());

    // use wisdomfile, if it exists
    if (fftw_import_wisdom_from_filename(wisdompath.c_str()) != 0) {
        plan = fftw_plan_dft_1d(n,in,out,sign,FFTW_WISDOM_ONLY|FFTW_PATIENT);
    }
    // if there was no wisdom (no or bad file), create some
    if (plan == nullptr) {
        plan = fftw_plan_dft_1d(n,in,out,sign,FFTW_PATIENT);
        fftw_export_wisdom_to_filename(wisdompath.c_str());
        vfps::Display::printText("Created some wisdom at "+wisdompath.str());
    }
    return plan;
}


fftwf_plan fft::prepareFFT(size_t n,
                                           fftwf_complex* in,
                                           fftwf_complex* out,
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

    std::stringstream wisdomfname;
    // find filename for wisdom
    wisdomfname << "wisdom_c" << dir << "c32_" << n << ".fftw";

    vfps::FSPath wisdompath(vfps::FSPath::datapath());
    wisdompath.append("fftwisdom/"+wisdomfname.str());

    // use wisdomfile, if it exists
    if (fftwf_import_wisdom_from_filename(wisdompath.c_str()) != 0) {
        plan = fftwf_plan_dft_1d(n,in,out,sign,FFTW_WISDOM_ONLY|FFTW_PATIENT);
    }
    // if there was no wisdom (no or bad file), create some
    if (plan == nullptr) {
        plan = fftwf_plan_dft_1d(n,in,out,sign,FFTW_PATIENT);
        fftwf_export_wisdom_to_filename(wisdompath.c_str());
        vfps::Display::printText("Created some wisdom at "+wisdompath.str());
    }
    return plan;
}

