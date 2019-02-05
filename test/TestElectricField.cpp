#include <boost/test/unit_test.hpp>

#include <vector>

#define INOVESA_ALLOW_PS_RESET 1

#include "defines.hpp"
#include "PS/ElectricField.hpp"
#include "PS/PhaseSpace.hpp"
#include "Z/ConstImpedance.hpp"



BOOST_AUTO_TEST_SUITE( ElectricField )

BOOST_AUTO_TEST_CASE( constructors ){
    constexpr uint32_t n = 32;

    vfps::PhaseSpace::resetSize(n,1);

    auto ps1 = std::make_shared<vfps::PhaseSpace>(
                -12,12,3,-12,12,4,nullptr,1,1);
    auto z1 = std::make_shared<vfps::ConstImpedance>(n, 1e9, 1);

    std::vector<uint32_t> buckets(1, 0);

    vfps::ElectricField f(ps1,z1, buckets, 0, nullptr, 1e6, 0.1,
                          1e-3, 1e9, 1e3, 1e-3);


    BOOST_CHECK_EQUAL_COLLECTIONS(f.getBuckets().data(),
                                  f.getBuckets().data()+f.getBuckets().size(),
                                  buckets.data(),
                                  buckets.data()+buckets.size());

    // TODO: add more checks
}

BOOST_AUTO_TEST_SUITE_END()
