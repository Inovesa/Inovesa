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
#include "Z/FreeSpaceCSR.hpp"


struct ElectricFieldFixture {
    explicit ElectricFieldFixture(std::vector<vfps::integral_t> filling,
                                  vfps::impedance_t Z)
      : spacing_bins(static_cast<size_t>(std::ceil(n*bunchspacing)))
      , spaced_bins(spacing_bins * filling.size())
      , filling(filling)
      , z(std::make_shared<vfps::ConstImpedance>(spaced_bins, 1e9, Z))
    {
        for (uint32_t i=0; i<filling.size(); i++) {
            if (filling[i] > 0) {
                buckets.push_back(static_cast<uint32_t>(filling.size())-1-i);
                bunches.emplace_back(filling[i]);
            }
        }

        vfps::PhaseSpace::resetSize(
                    n, static_cast<vfps::meshindex_t>(bunches.size()));

        ps = std::make_shared<vfps::PhaseSpace>(
                    -12,12,3,-12,12,4,nullptr,1,1,bunches);

        f = std::make_shared<vfps::ElectricField>(
                    ps, z, buckets, spacing_bins, nullptr,
                    1e6, 0.1, 1e-3, 1e9, 1e3, 1e-3);
    }

    static constexpr float bunchspacing = 8.0;

    static constexpr uint32_t n = 32;

    size_t spacing_bins;

    size_t spaced_bins;

    std::vector<uint32_t> buckets;

    std::vector<vfps::integral_t> bunches;

    std::vector<vfps::integral_t> filling;

    std::shared_ptr<vfps::PhaseSpace> ps;

    std::shared_ptr<vfps::Impedance> z;

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
    ElectricFieldFixture eff(T::pattern, 1);

    BOOST_CHECK_EQUAL_COLLECTIONS(eff.f->getBuckets().data(),
                                  eff.f->getBuckets().data()
                                  +eff.f->getBuckets().size(),
                                  eff.buckets.data(),
                                  eff.buckets.data()
                                  +eff.buckets.size());
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(padding , T, filling_patterns, T)
{
    ElectricFieldFixture eff(T::pattern, 1);

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
        for (vfps::meshindex_t x_bucket(0); x_bucket<eff.n;
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
 * @brief wake_potential test
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(wake_potential , T, filling_patterns, T)
{
    ElectricFieldFixture eff(T::pattern, 1);

    eff.f->padBunchProfiles();
    eff.f->wakePotential();

    auto padded_profile = eff.ps->getProjection(0);
    auto wake_potential = eff.f->getWakePotentials();
    auto scaling = eff.f->getWakeScaling();

    for (vfps::meshindex_t b=0; b<vfps::PhaseSpace::nb;b++) {
        for (vfps::meshindex_t x=0; x<eff.n; x++) {
            /* Shapes of profile and wake potential (DFT forth and back)
             * should match, as constant 1 is used as impedance.
             */
            BOOST_CHECK_SMALL(wake_potential[b][x]
                              -padded_profile[b][x]*scaling,
                              static_cast<vfps::integral_t>(1e-6));
        }
    }
}

/**
 * @brief forward_wake test whether FS CSR just influences forward wake
 *
 * This is actually ore of an integration test.
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(forward_wake, T, filling_patterns, T)
{
    ElectricFieldFixture eff(T::pattern, 0);

    eff.z->operator+=(vfps::FreeSpaceCSR(eff.spaced_bins, 1e6, 1e12));

    // original profiles
    eff.f->padBunchProfiles();
    eff.f->wakePotential();

    boost::multi_array<vfps::projection_t,2> profile_orig
            = eff.ps->getProjection(0);
    auto wakepot_orig = eff.f->getWakePotentials();

    // apply changes only in the second half
    constexpr vfps::meshindex_t change_from = eff.n/2;

    // modified profiles
    boost::multi_array<vfps::projection_t,2> profile_mod(profile_orig);
    for (vfps::meshindex_t b = 0; b < vfps::PhaseSpace::nb; b++) {
        auto delta = static_cast<vfps::projection_t>(0.1)
                     * profile_mod[b][change_from];
        profile_mod[b][change_from] -= delta;
        profile_mod[b][change_from+1] += delta;
        profile_mod[b][change_from+2] += delta;
        profile_mod[b][change_from+3] -= delta;
        eff.ps->setProjection(0,b,profile_mod[b]);
    }
    eff.f->padBunchProfiles();
    eff.f->wakePotential();
    auto wakepot_mod = eff.f->getWakePotentials();

    for (vfps::meshindex_t b = 0; b < vfps::PhaseSpace::nb; b++) {
        for (vfps::meshindex_t x = 0; x < change_from; x++) {
            BOOST_CHECK_SMALL(wakepot_orig[b][x]-wakepot_mod[b][x],
                              static_cast<vfps::integral_t>(1e-7));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
