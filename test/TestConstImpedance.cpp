// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#include <boost/test/unit_test.hpp>

#include<vector>
#include <iostream>

#include "defines.hpp"
#include "Z/ConstImpedance.hpp"

BOOST_AUTO_TEST_SUITE( ConstImpedance )

BOOST_AUTO_TEST_CASE( constructors ){
    constexpr size_t n = 6;
    vfps::impedance_t Z = 1e3;
    std::vector<vfps::impedance_t> zv1(n/2, Z);
    zv1.resize(n,0);

    vfps::ConstImpedance z1(n, 1e9, Z);

    BOOST_CHECK_EQUAL_COLLECTIONS(zv1.data(),zv1.data()+zv1.size(),
                                  z1.data(),z1.data()+n);
}

BOOST_AUTO_TEST_SUITE_END()
