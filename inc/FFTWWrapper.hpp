// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#pragma once

#include <algorithm>
#include <fftw3.h>
#include <memory>

#include "defines.hpp"

namespace fft {

typedef fftwf_plan plan;
typedef fftwf_complex complex;


enum class fft_direction : uint_fast8_t {
    forward, backward
};

complex* fft_alloc_complex(size_t n);
vfps::integral_t* fft_alloc_real(size_t n);

/**
 * @brief fft_cleanup to be implemented
 */
inline void fft_cleanup()
    {}

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


fftw_plan prepareFFT(size_t n, double* in, fftw_complex* out);

inline fftw_plan prepareFFT(size_t n, double* in, std::complex<double>* out)
    { return prepareFFT(n,in, reinterpret_cast<fftw_complex*>(out)); }

/**
 * @brief prepareFFT real to complex (Hermitian) FFT, "forward"
 * @param n
 * @param in
 * @param out
 * @return
 */
fftwf_plan prepareFFT(size_t n, float* in, fftwf_complex* out);

inline fftwf_plan prepareFFT(size_t n, float* in, std::complex<float>* out)
    { return prepareFFT(n,in, reinterpret_cast<fftwf_complex*>(out)); }

/**
 * @brief prepareFFT Hermitian (complex) to real FFT, "backward"
 * @param n
 * @param in
 * @param out
 * @return
 */
fftwf_plan prepareFFT(size_t n, fftwf_complex* in, float* out);

inline fftwf_plan prepareFFT(size_t n, std::complex<float>* in, float* out)
    {return prepareFFT(n,reinterpret_cast<fftwf_complex*>(in), out); }

/**
 * @brief prepareFFT (unmaintained for C2C FFT)
 */
fftw_plan prepareFFT(size_t n, fftw_complex *in,
                     fftw_complex *out,
                     fft_direction direction);

inline fftw_plan prepareFFT(size_t n, std::complex<double>* in,
                            std::complex<double>* out,
                            fft_direction direction)
    { return prepareFFT(n,reinterpret_cast<fftw_complex*>(in),
                        reinterpret_cast<fftw_complex*>(out),direction); }

/**
 * @brief prepareFFT (unmaintained for C2C FFT)
 */
fftwf_plan prepareFFT(size_t n, fftwf_complex* in,
                      fftwf_complex* out,
                      fft_direction direction);

inline fftwf_plan prepareFFT(size_t n, std::complex<float>* in,
                            std::complex<float>* out,
                            fft_direction direction)
    { return prepareFFT(n,reinterpret_cast<fftwf_complex*>(in),
                        reinterpret_cast<fftwf_complex*>(out),direction); }

} // namespace fftwwrapper
