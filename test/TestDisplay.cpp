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
    BOOST_CHECK_NO_THROW(vfps::make_display("foo.conf",false));
}

BOOST_AUTO_TEST_SUITE_END()
