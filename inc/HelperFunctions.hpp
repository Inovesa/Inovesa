// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include <ctime>
#include <iomanip>
#include <memory>
#include <mutex>
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
        const std::string& branch = GIT_BRANCH,
        const std::string& commit = GIT_COMMIT);

bool isOfFileType(std::string ending, std::string fname);

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

/**
 * @brief localtime threading-save wrapper for std::localtime
 * @param time time stamp to convert
 * @return struct holdig calendar date and time
 */
inline std::tm localtime(std::time_t time)
{
    std::tm bt {};
    #if defined(__STDC_LIB_EXT1__) // C++11 extension for safe alternatives
    localtime_s(&bt, &time);
    #elif defined(__unix__)
    localtime_r(&time, &bt); // unix "restartable" multi-threading alternative
    #else // fallback if no thread-safe alternative exists
    static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    bt = *std::localtime(&time);
    #endif
    return bt;
}

const std::string status_string(std::shared_ptr<PhaseSpace> ps,
                                float roatation,
                                float total_rotations);

uint64_t upper_power_of_two(uint64_t v);

} // namespace vfps

