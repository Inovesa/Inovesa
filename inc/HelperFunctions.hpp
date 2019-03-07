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

/** \file
 *  \brief Definitions of helper functions
 */

namespace vfps
{

const std::string copyright_notice() noexcept;

const std::string inovesa_version(
        const bool verbose=false,
        const int_fast16_t v_mayor = INOVESA_VERSION_MAJOR,
        const int_fast16_t v_minor = INOVESA_VERSION_MINOR,
        const int_fast16_t v_fix = INOVESA_VERSION_FIX,
        const std::string branch = GIT_BRANCH,
        const std::string commit = GIT_COMMIT);

bool isOfFileType(std::string ending, std::string fname);

const std::string status_string(std::shared_ptr<PhaseSpace> ps,
                                float roatation,
                                float total_rotations);

uint64_t upper_power_of_two(uint64_t v);

} // namespace vfps

