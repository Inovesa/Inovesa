// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#include <boost/test/unit_test.hpp>

#include <vector>

#define INOVESA_ALLOW_PS_RESET 1

#include "defines.hpp"
#include "PS/ElectricField.hpp"
#include "PS/PhaseSpace.hpp"
#include "Z/ConstImpedance.hpp"


struct ElectricFieldFixture {
    ElectricFieldFixture()
     : buckets(1, 0)
    {
        vfps::PhaseSpace::resetSize(n,1);

        ps = std::make_shared<vfps::PhaseSpace>(-12,12,3,-12,12,4,nullptr,1,1);
        z = std::make_shared<vfps::ConstImpedance>(n, 1e9, 1);

        f = std::make_shared<vfps::ElectricField>(
                    ps, z, buckets, 0, nullptr, 1e6, 0.1, 1e-3, 1e9, 1e3, 1e-3);
    }

    static constexpr uint32_t n = 32;

    std::vector<uint32_t> buckets;

    std::shared_ptr<vfps::PhaseSpace> ps;

    std::shared_ptr<vfps::ConstImpedance> z;

    std::shared_ptr<vfps::ElectricField> f;
};


constexpr uint32_t ElectricFieldFixture::n;

BOOST_FIXTURE_TEST_SUITE( ElectricField, ElectricFieldFixture )

BOOST_AUTO_TEST_CASE( constructors ){
    BOOST_CHECK_EQUAL_COLLECTIONS(f->getBuckets().data(),
                                  f->getBuckets().data()+f->getBuckets().size(),
                                  buckets.data(),
                                  buckets.data()+buckets.size());

    // TODO: add more checks
}

BOOST_AUTO_TEST_SUITE_END()
