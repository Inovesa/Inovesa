#define BOOST_TEST_MODULE messagestrings
#include <boost/test/unit_test.hpp>

#include <boost/algorithm/string.hpp>
#include <vector>


#include "MessageStrings.hpp"

void longLinesCheck(std::string text) {
    std::vector<std::string> lines;
    boost::split(lines, text, boost::is_any_of("\n"));
    for(auto line : lines) {
        // lines should be no longer than 80 characters
        BOOST_WARN_GT(line.size(), static_cast<std::size_t>(80));
    }
}

BOOST_AUTO_TEST_CASE(copyright_notice) {
    auto s = vfps::copyright_notice();

    // there should be at least five lines
    BOOST_CHECK_GE(std::count(s.begin(), s.end(), '\n'), static_cast<std::size_t>(5));

    longLinesCheck(s);
}

BOOST_AUTO_TEST_CASE(inovesa_version) {
    auto s = vfps::inovesa_version(true);
    longLinesCheck(s);
}
