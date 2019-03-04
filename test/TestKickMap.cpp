#include <boost/test/unit_test.hpp>

#include <array>
#include <boost/math/constants/constants.hpp>

#include "defines.hpp"
#define INOVESA_ALLOW_PS_RESET 1
#include "PS/PhaseSpace.hpp"
#include "SM/KickMap.hpp"

BOOST_AUTO_TEST_SUITE( KickMap )

BOOST_AUTO_TEST_CASE( linear_noclamp_map_x ){
    std::vector<vfps::integral_t> buckets{{1}};
    vfps::PhaseSpace::resetSize(4,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.0, 0.1, 0.0, 1.0,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.1, 0.7, 0.3, 0.0,
                                          0.9, 0.0, 0.0, 0.0 }};

    std::vector<vfps::meshaxis_t> offset{{ -1, 1.0, 0.0, 1.0 }};

    std::vector<vfps::meshdata_t> data2{{ 0.0, 0.3, 0.0, 0.0,
                                          0.0, 0.7, 0.7, 0.0,
                                          0.0, 0.0, 0.3, 0.0,
                                          0.1, 0.0, 0.0, 0.0 }};

    auto ps1 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1,data1.data());

    auto ps2 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1);

    vfps::KickMap km(
                ps1,ps2,4,4,buckets.size(),
                vfps::SourceMap::InterpolationType::linear,false,
                vfps::KickMap::Axis::x,nullptr);

    km.swapForce(offset);

    km.apply();

    for (auto i=0; i<16; i++) {
        BOOST_CHECK_SMALL(ps2->getData()[i]-data2[i],
                          static_cast<vfps::meshdata_t>(1e-5));
    }
}

BOOST_AUTO_TEST_CASE( linear_noclamp_map_y ){
    std::vector<vfps::integral_t> buckets{{1}};
    vfps::PhaseSpace::resetSize(4,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.0, 0.1, 0.0, 0.3,
                                          0.0, 0.3, 0.7, 0.0,
                                          0.1, 0.7, 0.3, 0.0,
                                          0.5, 0.0, 0.6, 0.0 }};

    std::vector<vfps::meshaxis_t> offset{{ -1, 1.0, 0.0, 1.0 }};

    std::vector<vfps::meshdata_t> data2{{ 0.0, 0.0, 0.1, 0.0,
                                          0.3, 0.7, 0.0, 0.0,
                                          0.1, 0.7, 0.3, 0.0,
                                          0.0, 0.6, 0.0, 0.0 }};

    auto ps1 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1,data1.data());

    auto ps2 = std::make_shared<vfps::PhaseSpace>(-1,1,2,-1,1,4,nullptr,1,1,
                                                  buckets,1);

    vfps::KickMap km(
                ps1,ps2,4,4,buckets.size(),
                vfps::SourceMap::InterpolationType::linear,false,
                vfps::KickMap::Axis::y,nullptr);

    km.swapForce(offset);

    km.apply();

    for (auto i=0; i<16; i++) {
        BOOST_CHECK_SMALL(ps2->getData()[i]-data2[i],
                          static_cast<vfps::meshdata_t>(1e-5));
    }
}

BOOST_AUTO_TEST_SUITE_END()
