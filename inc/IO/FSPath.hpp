/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
 *                                                                            *
 * This file is part of Inovesa.                                              *
 * Inovesa is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Inovesa is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           *
 ******************************************************************************/

#ifndef FSPATH_HPP
#define FSPATH_HPP

#include <string>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace vfps
{

/**
 * @brief The FSPath class wraps full paths to files and directories.
 */
class FSPath
{
public:
    FSPath(std::string path);

    /**
     * @brief append
     * @param path
     * @return *this
     */
    FSPath& append(std::string path);

    inline const std::string& str() const
        { return _path.string(); }

    inline const char* c_str() const
        { return _path.c_str(); }

    static std::string datapath();

private:
    static std::string expand_user(std::string path);

    void checkDirectory(std::string path);

    fs::path _path;
};

} // namespace vfps

#endif // FSPATH_HPP

