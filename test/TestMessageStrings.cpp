#define BOOST_TEST_MODULE messagestrings
#include <boost/test/unit_test.hpp>

#include <boost/algorithm/string.hpp>
#include <vector>


#include "MessageStrings.hpp"

BOOST_AUTO_TEST_CASE( copyright_notice ) {
    std::string s = vfps::copyright_notice();

    // there should be at least five lines
    BOOST_CHECK_GE(std::count(s.begin(), s.end(), '\n'), static_cast<std::size_t>(5));

    std::vector<std::string> lines;
    boost::split(lines, s, boost::is_any_of("\n"));
    for(auto line : lines) {
        // lines should be no longer than 80 characters
        BOOST_WARN_GT(line.size(), static_cast<std::size_t>(80));
    }
}
