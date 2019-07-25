// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include <string>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace vfps
{

/**
 * @brief wrapper for Windows compatibility of setenv
 *
 * @param name
 * @param value
 * @param overwrite
 */
int setenv(const char* name, const char* value, int overwrite);


/**
 * @brief wrapper for Windows compatibility of unsetenv
 *
 * @param name
 * @param value
 * @param overwrite
 */
int unsetenv(const char* name);

/**
 * @brief The FSPath class wraps full paths to files and directories.
 */
class FSPath
{
public:
    FSPath() = delete;

    explicit FSPath(std::string path);

    /**
     * @brief append
     * @param path
     * @return *this
     */
    FSPath& append(std::string path);

    inline const std::string& str() const
        { return _path.string(); }

    inline const char* c_str() const
        { return str().c_str(); }

    static std::string datapath();

    static std::string expand_user(std::string path);

    static bool validateDirectory(fs::path path);

private:
    fs::path _path;
};

} // namespace vfps

