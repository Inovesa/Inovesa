#include <boost/test/unit_test.hpp>

#include <array>
#include <boost/math/constants/constants.hpp>

#include "defines.hpp"
#define INOVESA_ALLOW_PS_RESET 1
#include "PS/PhaseSpace.hpp"
#include "SM/RotationMap.hpp"


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

    BOOST_CHECK_EQUAL(std::memcmp(data2.data(),ps2->getData(),data2.size()), 0);
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

    BOOST_CHECK_EQUAL(std::memcmp(data2.data(),ps2->getData(),data2.size()), 0);
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

    BOOST_CHECK_EQUAL(std::memcmp(data2.data(),ps2->getData(),data2.size()), 0);

    vfps::PhaseSpace::Position p0{0,0};
    auto p1 = rm.apply(p0);

    // expected result from the above operation
    vfps::PhaseSpace::Position p2{3,0};

    // chack whether difference is small (BOOST_CHECK_CLOSE checks relative)
    BOOST_CHECK_SMALL(p1.x-p2.x, static_cast<vfps::meshaxis_t>(1e-3f));
    BOOST_CHECK_SMALL(p1.y-p2.y, static_cast<vfps::meshaxis_t>(1e-3f));
}
