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
/**
 * @brief makePSFromPNG
 * @param fname file name to load phase space from
 * @param qmin (see PhaseSpace)
 * @param qmax (see PhaseSpace)
 * @param xscale (see PhaseSpace)
 * @param pmin (see PhaseSpace)
 * @param pmax (see PhaseSpace)
 * @param yscale (see PhaseSpace)
 * @param oclh (see PhaseSpace)
 * @param beam_charge (see PhaseSpace)
 * @param beam_current (see PhaseSpace)
 * @return pointer to loaded (and re-normalized) PhaseSpace
 *
 * @todo Add support for multi-bunch.
 */
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

std::unique_ptr<PhaseSpace> makePSFromTXT( const std::string& fname
                                         , int64_t ps_size
                                         , meshaxis_t qmin, meshaxis_t qmax
                                         , meshaxis_t pmin, meshaxis_t pmax
                                         , oclhptr_t oclh
                                         , const double beam_charge
                                         , const double beam_current
                                         , double xscale, double yscale
                                         );

/**
 * @brief saveToImage saves phase spaces to image
 * @param ps phasespace to be saved
 * @param ofname name to save data to (file ending defines type)
 *
 * Known file types are:
 * - PNG
 *
 * Data is saved 1:1 grid-point to pixel,
 * while the horizontal dimension gives space
 * and energy is aligned vertically.
 * Multi-bunch configurations will be saved side by side
 * in the space dimensions without any separation.
 */
void saveToImage(const PhaseSpace& ps,
                 const std::string& ofname);

} // namespace vfps

