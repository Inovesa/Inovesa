// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#include <stdlib.h>

#include "IO/FSPath.hpp"

vfps::FSPath::FSPath(std::string path)
  : _path(expand_user(path))
{
    validateDirectory(_path);
}

vfps::FSPath& vfps::FSPath::append(std::string path)
{
    _path /= path;
    validateDirectory(_path);
    return *this;
}

/* References:
 * https://stackoverflow.com/questions/4891006/how-to-create-a-folder-in-the-home-directory
 */
std::string vfps::FSPath::expand_user(std::string path)
{
  if (!path.empty() && path[0] == '~') {
    assert(path.size() == 1 || path[1] == '/');  // or other error handling
    char const* home = std::getenv("HOME");
    if (home || ((home = std::getenv("USERPROFILE")))) {
      path.replace(0, 1, home);
    }
    else {
      char const *hdrive = std::getenv("HOMEDRIVE");
      char const *hpath = std::getenv("HOMEPATH");
      assert(hdrive);  // or other error handling
      assert(hpath);
      path.replace(0, 1, std::string(hdrive) + hpath);
    }
  }
  return path;
}

bool vfps::FSPath::validateDirectory(boost::filesystem::path path)
{
    if (!fs::exists(path)) {
        if ((path.string().back()) == '/') {
            return fs::create_directories(path);
        } else {
            return fs::create_directories(path.parent_path());
        }
    }
    return false;
}

/* References:
 * https://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html#variables
 *
 * If XDG_DATA_HOME is either not set or empty,
 * a default equal to $HOME/.local/share should be used.
 */
std::string vfps::FSPath::datapath()
{
    char const *envstr = std::getenv("XDG_DATA_HOME");
    std::string rval;
    if (envstr != nullptr) {
        rval = std::string(envstr);
    } else {
        rval = expand_user("~/.local/share");
    }
    rval += "/inovesa/";
    return rval;
}

