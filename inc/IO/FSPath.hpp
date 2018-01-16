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

