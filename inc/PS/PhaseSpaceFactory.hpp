// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include <memory>
#include <string>

#include "PS/PhaseSpace.hpp"

namespace vfps {

#if INOVESA_USE_HDF5 == 1
std::unique_ptr<PhaseSpace> makePSFromHDF5( std::string fname
                                          , int64_t startdiststep
                                          , meshaxis_t qmin, meshaxis_t qmax
                                          , meshaxis_t pmin, meshaxis_t pmax
                                          , oclhptr_t oclh
                                          , const double beam_charge
                                          , const double beam_current
                                          , double xscale, double yscale
                                          );
#endif // INOVESA_USE_HDF5

std::unique_ptr<PhaseSpace> makePSFromPNG( std::string fname
                                         , meshaxis_t qmin, meshaxis_t qmax
                                         , meshaxis_t pmin, meshaxis_t pmax
                                         , oclhptr_t oclh
                                         , const double beam_charge
                                         , const double beam_current
                                         , double xscale, double yscale
                                         );

std::unique_ptr<PhaseSpace> makePSFromTXT( std::string fname, int64_t ps_size
                                         , meshaxis_t qmin, meshaxis_t qmax
                                         , meshaxis_t pmin, meshaxis_t pmax
                                         , oclhptr_t oclh
                                         , const double beam_charge
                                         , const double beam_current
                                         , double xscale, double yscale
                                         );


} // namespace vfps
