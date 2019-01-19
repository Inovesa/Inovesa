#include <boost/test/unit_test.hpp>

#include <array>
#include <boost/math/constants/constants.hpp>

#include "defines.hpp"
#define INOVESA_ALLOW_PS_RESET 1
#include "PS/PhaseSpace.hpp"
#include "SM/RotationMap.hpp"


BOOST_AUTO_TEST_CASE( centered90 ){
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
                vfps::RotationMap::cubic,false,0,nullptr);

    rm.apply();

    BOOST_CHECK_EQUAL(std::memcmp(data2.data(),ps2->getData(),data2.size()), 0);
}
