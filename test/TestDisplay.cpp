// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Fistname Lastename
 */

#include <boost/test/unit_test.hpp>

#include <string>

#include "defines.hpp"
#include "IO/Display.hpp"

BOOST_AUTO_TEST_SUITE( Display )

const std::string errormessage("some error message");

bool checkMessage(const std::exception& ex)
{
    BOOST_CHECK_EQUAL(ex.what(), errormessage);
    return true;
}

[[noreturn]] void throwDisplayException() {
    throw vfps::DisplayException(errormessage);
}

BOOST_AUTO_TEST_CASE( DisplayException ){
    BOOST_CHECK_EXCEPTION(throwDisplayException(),
                          vfps::DisplayException,
                          checkMessage);
}

BOOST_AUTO_TEST_CASE( MakeDisplay ) {
    vfps::Display::silent_mode = true;
    BOOST_CHECK_NO_THROW(vfps::make_display("foo.conf",false));
    vfps::Display::silent_mode = false;
}

BOOST_AUTO_TEST_CASE( AbortHandling ) {
    // abourt shoud be initialized with false
    constexpr auto original_state = false;
    BOOST_CHECK(vfps::Display::abort == original_state);

    // SIGINT_handler sets abort to true
    BOOST_CHECK_NO_THROW(vfps::Display::SIGINT_handler(SIGINT));
    BOOST_CHECK(vfps::Display::abort);

    // calling multiple times does not change this again
    BOOST_CHECK_NO_THROW(vfps::Display::SIGINT_handler(SIGINT));
    BOOST_CHECK(vfps::Display::abort);

    // restore original value (falce)
    vfps::Display::abort = original_state;
}

BOOST_AUTO_TEST_SUITE_END()
