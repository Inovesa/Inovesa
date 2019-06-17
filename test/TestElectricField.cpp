// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/vector.hpp>
#include <vector>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::one_div_root_two_pi;

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

        spacing_bins = static_cast<size_t>(std::ceil(n*bunchspacing));

        spaced_bins = spacing_bins * filling.size();

        vfps::PhaseSpace::resetSize(n, bunches.size());

        ps = std::make_shared<vfps::PhaseSpace>(
                    -12,12,3,-12,12,4,nullptr,1,1,bunches);

        // impedance with constant value of 1
        z = std::make_shared<vfps::ConstImpedance>(spaced_bins, 1e9, 1);

        f = std::make_shared<vfps::ElectricField>(
                    ps, z, buckets, spacing_bins, nullptr,
                    1e6, 0.1, 1e-3, 1e9, 1e3, 1e-3);
    }

    static constexpr uint32_t n = 32;

    size_t spaced_bins;

    size_t spacing_bins;

    std::vector<uint32_t> buckets;

    std::vector<vfps::integral_t> bunches;

    std::vector<vfps::integral_t> filling;

    std::shared_ptr<vfps::PhaseSpace> ps;

    std::shared_ptr<vfps::ConstImpedance> z;

    std::shared_ptr<vfps::ElectricField> f;

    static constexpr vfps::integral_t zoom = 1.0;
    static constexpr vfps::integral_t zoom2 = zoom*zoom;
};

constexpr uint32_t ElectricFieldFixture::n;
constexpr vfps::integral_t ElectricFieldFixture::zoom;
constexpr vfps::integral_t ElectricFieldFixture::zoom2;

struct filling {
    explicit filling(std::vector<vfps::integral_t> pattern)
      : pattern(pattern) {}

    std::vector<vfps::integral_t> pattern;
};

struct filling1 : public filling {
    filling1() : filling({1.0}) {}
};

struct filling2 : public filling {
    filling2() : filling({0.75, 0, 0.25}) {}
};

typedef boost::mpl::vector<filling1,
                           filling2
                          > filling_patterns;


BOOST_AUTO_TEST_SUITE( ElectricField )

BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructors , T, filling_patterns, T)
{
    ElectricFieldFixture eff(T::pattern);

    BOOST_CHECK_EQUAL_COLLECTIONS(eff.f->getBuckets().data(),
                                  eff.f->getBuckets().data()
                                  +eff.f->getBuckets().size(),
                                  eff.buckets.data(),
                                  eff.buckets.data()
                                  +eff.buckets.size());
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(padding , T, filling_patterns, T)
{
    ElectricFieldFixture eff(T::pattern);

    eff.ps->updateXProjection();
    eff.f->padBunchProfiles();

    std::vector<vfps::integral_t> correct_solution(eff.spaced_bins,0);
    size_t bucket(0);
    vfps::meshindex_t x(0);
    auto axis = eff.ps->getAxis(0);
    /* Time and space coordinates are anti-parallel.
     * As pattern is in time but grid is in space,
     * we have to interate reversed.
     */
    for(auto norm = T::pattern.rbegin(); norm != T::pattern.rend(); ++norm)
    {
        x = bucket*eff.spacing_bins;
        for (vfps::meshindex_t x_bucket=0; x_bucket<eff.n;
             x_bucket++, x++) {
            correct_solution[x]
                  =  *norm * one_div_root_two_pi<vfps::integral_t>()
                  * std::exp((-0.5)*axis->at(x_bucket)
                             * axis->at(x_bucket)/eff.zoom2);
        }
        bucket++;
    }

    auto test_solution = eff.f->getPaddedBunchProfiles();
    for (size_t x=0; x<eff.spaced_bins; x++) {
        BOOST_CHECK_CLOSE(correct_solution[x], test_solution[x], 1e-4);
    }
}


/**
 * @brief WakePotential test
 *
 * @todo check normalization
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(WakePotential , T, filling_patterns, T)
{
    ElectricFieldFixture eff(T::pattern);

    eff.ps->updateXProjection();
    eff.f->padBunchProfiles();
    eff.f->wakePotential();

    auto padded_profile = eff.ps->getProjection(0);
    auto wake_potential = eff.f->getWakePotentials();

    for (vfps::meshindex_t b=0; b<vfps::PhaseSpace::nb;b++) {
        for (vfps::meshindex_t x=0; x<eff.n; x++) {
            /* Shapes of profile and wake potential (DFT forth and back)
             * should match, as constant 1 is used as impedance.
             */
            BOOST_CHECK_SMALL(wake_potential[b][x]/wake_potential[b][eff.n/2]
                              -padded_profile[b][x]/padded_profile[b][eff.n/2],
                              static_cast<vfps::integral_t>(1e-6));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
