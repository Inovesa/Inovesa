// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/vector.hpp>

#include <vector>

#define INOVESA_ALLOW_PS_RESET 1

#include "defines.hpp"
#include "PS/ElectricField.hpp"
#include "PS/PhaseSpace.hpp"
#include "Z/ConstImpedance.hpp"


struct ElectricFieldFixture {
    explicit ElectricFieldFixture(std::vector<vfps::integral_t> filling)
      : filling(filling)
    {
        for (uint32_t i=0; i<filling.size(); i++) {
            if (filling[i] > 0) {
                buckets.push_back(static_cast<uint32_t>(filling.size())-1-i);
                bunches.emplace_back(filling[i]);
            }
        }
        constexpr double bunchspacing = 4.0;

        size_t spaced_bins = static_cast<size_t>(
                    std::ceil(n*buckets.size()*bunchspacing));

        vfps::PhaseSpace::resetSize(n,1);

        ps = std::make_shared<vfps::PhaseSpace>(
                    -12,12,3,-12,12,4,nullptr,1,1,bunches);
        z = std::make_shared<vfps::ConstImpedance>(spaced_bins, 1e9, 1);

        f = std::make_shared<vfps::ElectricField>(
                    ps, z, buckets, 0, nullptr,
                    1e6, 0.1, 1e-3, 1e9, 1e3, 1e-3);
    }

    static constexpr uint32_t n = 32;

    std::vector<uint32_t> buckets;

    std::vector<vfps::integral_t> bunches;

    std::vector<vfps::integral_t> filling;

    std::shared_ptr<vfps::PhaseSpace> ps;

    std::shared_ptr<vfps::ConstImpedance> z;

    std::shared_ptr<vfps::ElectricField> f;
};

constexpr uint32_t ElectricFieldFixture::n;

struct ElectricFieldFixture1 : public ElectricFieldFixture{
    ElectricFieldFixture1()
      : ElectricFieldFixture(std::vector<vfps::integral_t>{1.0})
    {}
};

struct ElectricFieldFixture2 : public ElectricFieldFixture{
    ElectricFieldFixture2()
      : ElectricFieldFixture(std::vector<vfps::integral_t>{0.75, 0, 0.25})
    {}
};

typedef boost::mpl::vector<ElectricFieldFixture1,
                           ElectricFieldFixture2
                          > ElectricFieldFixtures;


BOOST_AUTO_TEST_SUITE( ElectricField )

BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructors , T, ElectricFieldFixtures, T)
{
    BOOST_CHECK_EQUAL_COLLECTIONS(T::f->getBuckets().data(),
                                  T::f->getBuckets().data()
                                  +T::f->getBuckets().size(),
                                  T::buckets.data(),
                                  T::buckets.data()
                                  +T::buckets.size());
}

BOOST_AUTO_TEST_SUITE_END()
