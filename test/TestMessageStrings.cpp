#include <boost/test/unit_test.hpp>

#include <boost/algorithm/string.hpp>
#include <vector>


#include "MessageStrings.hpp"

// normaly, there should be no long lines
void longLinesCheck(std::string text) {
    std::vector<std::string> lines;
    boost::split(lines, text, boost::is_any_of("\n"));
    for(auto line : lines) {
        // lines should be no longer than 80 characters
        BOOST_CHECK_LE(line.size(), static_cast<std::size_t>(80));
    }
}

BOOST_AUTO_TEST_CASE(copyright_notice) {
    auto s = vfps::copyright_notice();

    // there should be at least five lines
    BOOST_CHECK_GE(std::count(s.begin(), s.end(), '\n'), static_cast<std::size_t>(5));

    longLinesCheck(s);
}

BOOST_AUTO_TEST_CASE(inovesa_version) {
    longLinesCheck(vfps::inovesa_version(true,1,2,3,"master","$(hash)"));

    BOOST_CHECK_EQUAL(vfps::inovesa_version(false,1,2,3,"master","$(hash)"),
                      "v1.2.3");
    BOOST_CHECK_EQUAL(vfps::inovesa_version(false,1,3,0,"v1.3","$(hash)"),
                      "v1.3.0, Commit: $(hash)");
    BOOST_CHECK_EQUAL(vfps::inovesa_version(false,1,3,0,"$(feature)","$(hash)"),
                      "Branch: $(feature), Commit: $(hash)");
    BOOST_CHECK_EQUAL(vfps::inovesa_version(false,1,3,-1,"v1.3","$(hash)"),
                      "v1.3 alpha, Commit: $(hash)");
    BOOST_CHECK_EQUAL(vfps::inovesa_version(false,1,3,-2,"v1.3","$(hash)"),
                      "v1.3 beta, Commit: $(hash)");
    BOOST_CHECK_EQUAL(vfps::inovesa_version(false,1,3,-3,"v1.3","$(hash)"),
                      "v1.3 RC1, Commit: $(hash)");
}

BOOST_AUTO_TEST_CASE(status_string) {
    vfps::PhaseSpace::setSize(32,1);
    auto ps = std::make_shared<vfps::PhaseSpace>(-12,12,1,-32,32,2,nullptr,1,1);

    (*ps)[0][16][16] = 10;
    ps->updateXProjection();
    ps->integrate();
    auto s = vfps::status_string(ps,0,1);
    longLinesCheck(s);
    BOOST_CHECK_EQUAL(s[17],'-');

    (*ps)[0][16][16] = -10;
    ps->updateXProjection();
    ps->integrate();
    s = vfps::status_string(ps,0,1);
    longLinesCheck(s);
    BOOST_CHECK_EQUAL(s[17],'+');
}
