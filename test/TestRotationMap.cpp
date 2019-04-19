// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 */

#include <boost/test/unit_test.hpp>

#include <array>
#include <boost/math/constants/constants.hpp>

#include "defines.hpp"
#define INOVESA_ALLOW_PS_RESET 1
#include "PS/PhaseSpace.hpp"
#include "SM/RotationMap.hpp"

BOOST_AUTO_TEST_SUITE( RotationMap )

BOOST_AUTO_TEST_CASE( linear_noclamp_map ){
    std::vector<vfps::integral_t> buckets{{1}};
    vfps::PhaseSpace::resetSize(4,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    std::vector<vfps::meshdata_t> data2{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    auto ps1 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1,data1.data());

    auto ps2 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1);

    vfps::RotationMap rm(
                ps1,ps2,4,4,boost::math::constants::half_pi<vfps::meshaxis_t>(),
                vfps::SourceMap::InterpolationType::linear,false,0,nullptr);

    rm.apply();

    for (auto i=0; i<16; i++) {
        BOOST_CHECK_SMALL(ps2->getData()[i]-data2[i],
                          static_cast<vfps::meshdata_t>(1e-5));
    }
}

BOOST_AUTO_TEST_CASE( linear_offcenter ){
    std::vector<vfps::integral_t> buckets{{1}};
    vfps::PhaseSpace::resetSize(4,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.1, 0.2, 0.0,
                                          0.0, 0.3, 0.4, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    std::vector<vfps::meshdata_t> data2{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.4, 0.3, 0.0,
                                          0.0, 0.2, 0.1, 0.0 }};

    auto ps1 = std::make_shared<vfps::PhaseSpace>(-2,1,2,-2,2,4,nullptr,1,1,
                                                  buckets,1,data1.data());

    auto ps2 = std::make_shared<vfps::PhaseSpace>(-2,1,2,-2,2,4,nullptr,1,1,
                                                  buckets,1);

    vfps::RotationMap rm(
                ps1,ps2,4,4,boost::math::constants::pi<vfps::meshaxis_t>(),
                vfps::SourceMap::InterpolationType::linear,false,0,nullptr);

    rm.apply();

    for (auto i=0; i<16; i++) {
        BOOST_CHECK_SMALL(ps2->getData()[i]-data2[i],
                          static_cast<vfps::meshdata_t>(1e-5));
    }
}

BOOST_AUTO_TEST_CASE( quadratic_clamp_nomap ){
    std::vector<vfps::integral_t> buckets{{1}};
    vfps::PhaseSpace::resetSize(4,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    std::vector<vfps::meshdata_t> data2{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    auto ps1 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1,data1.data());

    auto ps2 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1);


    BOOST_CHECK_THROW(
            vfps::RotationMap(
                    ps1,ps2,4,4,
                    boost::math::constants::half_pi<vfps::meshaxis_t>(),
                    vfps::SourceMap::InterpolationType::quadratic,true,0,
                    nullptr),
            std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( quadratic_noclamp_nomap ){
    std::vector<vfps::integral_t> buckets{{1}};
    vfps::PhaseSpace::resetSize(4,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    std::vector<vfps::meshdata_t> data2{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    auto ps1 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1,data1.data());

    auto ps2 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1);

    vfps::RotationMap rm(
                ps1,ps2,4,4,boost::math::constants::half_pi<vfps::meshaxis_t>(),
                vfps::SourceMap::InterpolationType::quadratic,false,0,nullptr);

    rm.apply();

    for (auto i=0; i<16; i++) {
        BOOST_CHECK_SMALL(ps2->getData()[i]-data2[i],
                          static_cast<vfps::meshdata_t>(1e-5));
    }
}

BOOST_AUTO_TEST_CASE( cubic_clamp_nomap ){
    std::vector<vfps::integral_t> buckets{{1}};
    vfps::PhaseSpace::resetSize(4,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    std::vector<vfps::meshdata_t> data2{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    auto ps1 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1,data1.data());

    auto ps2 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1);
    BOOST_CHECK_THROW(
            vfps::RotationMap(
                    ps1, ps2, 4, 4,
                    boost::math::constants::half_pi<vfps::meshaxis_t>(),
                    vfps::SourceMap::InterpolationType::cubic, true,
                    0, nullptr),
            std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( cubic_clamp_map ){
    std::vector<vfps::integral_t> buckets{{1}};
    vfps::PhaseSpace::resetSize(4,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    std::vector<vfps::meshdata_t> data2{{ 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.7, 0.3, 0.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.0, 0.0, 0.0, 0.0 }};

    auto ps1 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1,data1.data());

    auto ps2 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1);

    vfps::RotationMap rm(
                ps1, ps2, 4, 4,
                boost::math::constants::half_pi<vfps::meshaxis_t>(),
                vfps::SourceMap::InterpolationType::cubic, true,
                data1.size(), nullptr);

    rm.apply();


    for (auto i=0; i<16; i++) {
        BOOST_CHECK_SMALL(ps2->getData()[i]-data2[i],
                          static_cast<vfps::meshdata_t>(1e-5));
    }

    vfps::PhaseSpace::Position p1{0,0};
    rm.applyTo(p1);

    // expected result from the above operation
    vfps::PhaseSpace::Position p2{3,0};

    // chack whether difference is small (BOOST_CHECK_CLOSE checks relative)
    BOOST_CHECK_SMALL(p1.x-p2.x, static_cast<vfps::meshaxis_t>(1e-3f));
    BOOST_CHECK_SMALL(p1.y-p2.y, static_cast<vfps::meshaxis_t>(1e-3f));
}

BOOST_AUTO_TEST_SUITE_END()
