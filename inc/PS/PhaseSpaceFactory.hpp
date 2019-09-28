// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include <memory>
#include <string>

#include "PS/PhaseSpace.hpp"

namespace vfps {

#if INOVESA_USE_HDF5 == 1
std::unique_ptr<PhaseSpace> makePSFromHDF5( const std::string& fname
                                          , int64_t startdiststep
                                          , meshaxis_t qmin, meshaxis_t qmax
                                          , meshaxis_t pmin, meshaxis_t pmax
                                          , oclhptr_t oclh
                                          , const double beam_charge
                                          , const double beam_current
                                          , double xscale, double yscale
                                          );
#endif // INOVESA_USE_HDF5

#if INOVESA_USE_PNG == 1
std::unique_ptr<PhaseSpace> makePSFromPNG(const std::string& fname
                                         , meshaxis_t qmin, meshaxis_t qmax
                                         , double xscale
                                         , meshaxis_t pmin, meshaxis_t pmax
                                         , double yscale
                                         , oclhptr_t oclh
                                         , const double beam_charge
                                         , const double beam_current
                                         );
#endif // INOVESA_USE_PNG


void saveToImage( const PhaseSpace& ps,
                  const std::string& ofname);

std::unique_ptr<PhaseSpace> makePSFromTXT( const std::string fname, int64_t ps_size
                                         , meshaxis_t qmin, meshaxis_t qmax
                                         , meshaxis_t pmin, meshaxis_t pmax
                                         , oclhptr_t oclh
                                         , const double beam_charge
                                         , const double beam_current
                                         , double xscale, double yscale
                                         );

} // namespace vfps

