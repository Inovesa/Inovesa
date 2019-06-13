// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#include <boost/test/unit_test.hpp>

#include <vector>
#include <iostream>

#include "defines.hpp"
#include "Z/Impedance.hpp"

BOOST_AUTO_TEST_SUITE( Impedance )

BOOST_AUTO_TEST_CASE( constructors ){
    vfps::impedance_t j{0,0.5};

    std::vector<vfps::impedance_t> zv(8);
    zv[5] = j;

    vfps::Impedance z1(zv, 1e9);

    BOOST_CHECK_EQUAL(z1[0],static_cast<vfps::impedance_t>(0));
    BOOST_CHECK_EQUAL(z1[5],j);
    BOOST_CHECK_CLOSE(z1.getRuler()->max(), 1e9, 1e-3);

    vfps::Impedance z2(z1);

    BOOST_CHECK_EQUAL(z2[0],static_cast<vfps::impedance_t>(0));
    BOOST_CHECK_EQUAL(z2[5],j);
}

BOOST_AUTO_TEST_CASE( file_read_in ){
    constexpr int n = 6;
    std::vector<vfps::impedance_t> zv1(n);
    for (int32_t i = 0; i < n; i++) {
        zv1[i] = vfps::impedance_t(i-2, 2.5-i);
    }

    vfps::Impedance z1("impedance.dat", 1e9);

    BOOST_CHECK_EQUAL(z1.nFreqs(), n);
    BOOST_CHECK_EQUAL_COLLECTIONS(zv1.data(),zv1.data()+zv1.size(),
                                  z1.data(),z1.data()+n);
}

BOOST_AUTO_TEST_CASE( addition ){
    constexpr int n = 8;
    std::vector<vfps::impedance_t> zv1(n);
    std::vector<vfps::impedance_t> zv2(n);

    for (size_t i = 0; i < n; i++) {
        zv1[i] = vfps::impedance_t(i, -i+4);
        zv2[i] = vfps::impedance_t(i-2, 3);
    }

    vfps::Impedance z1(zv1, 1e9);
    vfps::Impedance z2(zv2, 1e9);
    vfps::Impedance z3(n,1e9);
    z3 = z1 + z2;

    BOOST_CHECK_EQUAL(z1.nFreqs(), n);
    BOOST_CHECK_EQUAL(z2.nFreqs(), n);
    BOOST_CHECK_EQUAL(z3.nFreqs(), n);

    for (size_t i = 0; i < n; i++) {
        BOOST_CHECK_EQUAL(z3[i],zv1[i]+zv2[i]);
    }

    z1 += z2;
    // z2 still unchanged, but z1 == z3
    for (size_t i = 0; i < n; i++) {
        BOOST_CHECK_EQUAL(z2[i],zv2[i]);
        BOOST_CHECK_EQUAL(z1[i],z3[i]);
    }

    vfps::Impedance z4(z1);
    z1.swap(z2);
    for (size_t i = 0; i < n; i++) {
        BOOST_CHECK_EQUAL(z1[i],zv2[i]);
        BOOST_CHECK_EQUAL(z2[i],z3[i]);
        BOOST_CHECK_EQUAL(z4[i],z2[i]);
    }
}

BOOST_AUTO_TEST_SUITE_END()
