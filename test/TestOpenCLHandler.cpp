// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <strstream>

#include "CL/OpenCLHandler.hpp"

BOOST_AUTO_TEST_SUITE( OpenCLHandler )
#if INOVESA_USE_OPENCL == 1

BOOST_AUTO_TEST_CASE( list_devices ){
    auto stdbuffer = std::cout.rdbuf();
    std::stringstream sstream;
    std::cout.rdbuf(sstream.rdbuf());
    BOOST_CHECK_NO_THROW(OCLH::listCLDevices());
    std::cout.rdbuf(stdbuffer);
}

#endif // INOVESA_USE_OPENCL
BOOST_AUTO_TEST_SUITE_END()
