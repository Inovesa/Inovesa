// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#include <boost/test/unit_test.hpp>

#include "defines.hpp"
#include "PS/Ruler.hpp"

BOOST_AUTO_TEST_SUITE( Ruler )

BOOST_AUTO_TEST_CASE( constructors ){
    vfps::Ruler<int32_t> axis(15,-4,4,{{"foo", 2.0}});
    BOOST_CHECK_EQUAL(axis.zerobin(), 7);

    // min > max
    BOOST_CHECK_THROW( vfps::Ruler<int32_t>(7,3,-3),
                       std::invalid_argument );

    auto axis2(axis);
    BOOST_CHECK_EQUAL(axis2.zerobin(), 7);
    BOOST_CHECK_EQUAL(axis2.length(), 8);
    BOOST_CHECK_EQUAL(axis2.scale("foo"), 2.0);
}

BOOST_AUTO_TEST_SUITE_END()
