// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include <iomanip>
#include <memory>
#include <sstream>


#include "PS/PhaseSpace.hpp"

namespace vfps
{

const std::string copyright_notice() noexcept;

/**
 * @brief inovesa_version
 * @param verbose
 * @return
 *
 * For branches leading to a release and for releases,
 * external applications rely on the format of the output
 * to determine the Inovesa feature level.
 * So, the string will always begin with "v{major}.{minor}"
 * followed by either ".{fix}" for releases
 * or " {descriptor}" for pre-releas versions.
 * Development versions do not have to follow this convention.
 */
const std::string inovesa_version(
        const bool verbose=false,
        const int_fast16_t v_mayor = INOVESA_VERSION_MAJOR,
        const int_fast16_t v_minor = INOVESA_VERSION_MINOR,
        const int_fast16_t v_fix = INOVESA_VERSION_FIX,
        const std::__cxx11::string branch = GIT_BRANCH,
        const std::__cxx11::string commit = GIT_COMMIT);

const std::string status_string(std::shared_ptr<PhaseSpace> ps,
                                float roatation,
                                float total_rotations);

} // namespace vfps

