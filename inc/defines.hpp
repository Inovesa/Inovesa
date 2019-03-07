// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include <complex>
#include <string>
#include <fftw3.h>

#include "InovesaConfig.hpp"

//#define INOVESA_SYNC_CL


/** \file
 *  \brief global definitions
 */

/** \namespace vfps
 *  \brief namespace Vlasov-Fokker-Planck solver
 */
namespace vfps
{
// has to be uint32_t (same as cl_uint) for OpenCL support
typedef uint32_t meshindex_t;

/* currently all of the below types have to be the same
 * use this to switch them all */
typedef float data_t;

typedef data_t csrpower_t;
typedef std::complex<csrpower_t> impedance_t;

typedef data_t frequency_t;
typedef data_t meshaxis_t;
typedef data_t meshdata_t;
typedef data_t interpol_t;
typedef data_t integral_t;
typedef data_t timeaxis_t;

typedef integral_t projection_t;

namespace physcons
{
/// speed of light (in m/s)
constexpr double c=2.99792458e8;

/// Vacuum permittivity (in F/m)
constexpr double epsilon0 = 8.854187817e-12;

/// Impedance of free space
constexpr double Z0 = 1/(epsilon0*c);

/// charge of one electron (in C)
constexpr double e=1.602e-19;

/// Alfven current (in A)
constexpr double IAlfven = 17045;

/// electron rest energy (in eV)
constexpr double me = 510998.9;

/// magnetic constant (in H/m)
constexpr double mu0 = 1/(epsilon0*c*c);
} // namespace physcons

} // namespace vfps

#ifndef INOVESA_USE_OPENGL
#define INOVESA_USE_OPENGL 0
#endif
#ifndef DEBUG
#define DEBUG 0
#endif
#ifndef INOVESA_ENABLE_INTERRUPT
#define INOVESA_ENABLE_INTERRUPT 0
#endif
#ifndef INOVESA_ENABLE_CLPROFILING
#define INOVESA_ENABLE_CLPROFILING 0
#endif
#ifndef INOVESA_USE_OPENCL
#define INOVESA_USE_OPENCL 0
#endif
#ifndef INOVESA_USE_CLFFT
#define INOVESA_USE_CLFFT 0
#endif
#ifndef INOVESA_USE_OPENGL
#define INOVESA_USE_OPENGL 0
#endif
#ifndef INOVESA_USE_HDF5
#define INOVESA_USE_HDF5 0
#endif
#ifndef INOVESA_USE_PNG
#define INOVESA_USE_PNG 0
#endif

