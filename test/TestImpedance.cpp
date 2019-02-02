#include <boost/test/unit_test.hpp>

#include<vector>
#include <iostream>

#include "defines.hpp"
#include "Z/Impedance.hpp"

BOOST_AUTO_TEST_CASE( impedance_constructors ){
    vfps::impedance_t j{0,1};

    std::vector<vfps::impedance_t> zv(8);
    zv[5] = j;

    vfps::Impedance Z1(zv, 1e9);

    BOOST_CHECK_EQUAL(Z1[0],static_cast<vfps::impedance_t>(0));
    BOOST_CHECK_EQUAL(Z1[5],j);

    vfps::Impedance Z2(Z1);

    BOOST_CHECK_EQUAL(Z2[0],static_cast<vfps::impedance_t>(0));
    BOOST_CHECK_EQUAL(Z2[5],j);
}

BOOST_AUTO_TEST_CASE( impedance_addition ){
    constexpr size_t n = 8;
    std::vector<vfps::impedance_t> zv1(n);
    std::vector<vfps::impedance_t> zv2(n);

    for (int i = 0; i < n; i++) {
        zv1[i] = vfps::impedance_t(i, -i+4);
        zv2[i] = vfps::impedance_t(i-2, 3);
    }

    vfps::Impedance z1(zv1, 1e9);
    vfps::Impedance z2(zv2, 1e9);
    vfps::Impedance z3 = z1 + z2;

    for (int i = 0; i < n; i++) {
        BOOST_CHECK_EQUAL(z3[i],zv1[i]+zv2[i]);
    }

    z1 += z2;
    // z2 still unchanged, but z1 == z3
    for (int i = 0; i < n; i++) {
        BOOST_CHECK_EQUAL(z2[i],zv2[i]);
        BOOST_CHECK_EQUAL(z1[i],z3[i]);
    }
}
