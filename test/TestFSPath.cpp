#include <boost/test/unit_test.hpp>

#include <vector>
#include <stdlib.h>

#include "defines.hpp"
#include "IO/FSPath.hpp"


BOOST_AUTO_TEST_SUITE( FSPath )

BOOST_AUTO_TEST_CASE( home_finder ){
    // save environment variables
    char const * homepath_ptr = std::getenv("HOME");
    const std::string homepath_orig(homepath_ptr ? homepath_ptr : "");
    char const * userprof_prt = std::getenv("USERPROFILE");
    const std::string userprof_orig(userprof_prt ? userprof_prt : "");

    // find using defaults
    std::string homepath;
    BOOST_CHECK_NO_THROW(homepath = vfps::FSPath::expand_user("~"));

    // home should always exist, so it is not created
    BOOST_CHECK_EQUAL(vfps::FSPath::validateDirectory(homepath), false);


    // find in HOME
    setenv("HOME",homepath.c_str(), 0);
    unsetenv("USERPROFILE");
    BOOST_CHECK_NO_THROW(homepath = vfps::FSPath::expand_user("~"));
    BOOST_CHECK_EQUAL(vfps::FSPath::validateDirectory(homepath), false);

    // find in USERPROFILE
    setenv("USERPROFILE",homepath.c_str(), 0);
    unsetenv("HOME");
    BOOST_CHECK_NO_THROW(homepath = vfps::FSPath::expand_user("~"));
    BOOST_CHECK_EQUAL(vfps::FSPath::validateDirectory(homepath), false);

    // restore environment variables
    setenv("HOME",homepath_orig.c_str(), 1);
    setenv("USERPROFILE",userprof_orig.c_str(), 1);
}

BOOST_AUTO_TEST_CASE( datapath ){
    std::string datapath;
    BOOST_CHECK_NO_THROW(datapath = vfps::FSPath::datapath());

    // datapath ends on "/inovesa/"
    std::string cmpstring = "/inovesa/";
    BOOST_CHECK_EQUAL_COLLECTIONS(
                datapath.data()+datapath.size()-cmpstring.size(),
                datapath.data()+datapath.size(),
                cmpstring.data(),
                cmpstring.data()+cmpstring.size());
}

BOOST_AUTO_TEST_CASE( constructor ){
    BOOST_CHECK_NO_THROW(vfps::FSPath("~"));
}

BOOST_AUTO_TEST_CASE( append_directory ){
    std::string randstring("UpTwf1eZhJE0cVNV");

    vfps::FSPath path("~");
    BOOST_CHECK_NO_THROW(path.append(randstring));

    auto pathstr = path.str();
    BOOST_CHECK_EQUAL_COLLECTIONS(
                pathstr.data()+pathstr.size()-randstring.size(),
                pathstr.data()+pathstr.size(),
                randstring.data(),
                randstring.data()+randstring.size());
}

BOOST_AUTO_TEST_SUITE_END()
