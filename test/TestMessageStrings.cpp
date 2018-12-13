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

BOOST_AUTO_TEST_CASE(status_string) {
    vfps::PhaseSpace::setSize(32,32,1);
    auto ps = std::make_shared<vfps::PhaseSpace>(32,-1,1,3,-2,2,4,nullptr, 1,1);

    auto s = vfps::status_string(ps,0,1);
    longLinesCheck(s);
}
