#include <boost/test/unit_test.hpp>

#include <array>

#include "defines.hpp"
#define INOVESA_ALLOW_PS_RESET 1
#include "PS/PhaseSpace.hpp"

BOOST_AUTO_TEST_CASE( phasespace_constructors ){
    {
    vfps::PhaseSpace::resetSize(32,1);
    // default constructor
    vfps::PhaseSpace ps1(-12,12,3,-12,12,4,nullptr,1,1);
    BOOST_CHECK_EQUAL(ps1.getMax(0),  12);
    BOOST_CHECK_EQUAL(ps1.getMin(0), -12);
    BOOST_CHECK_EQUAL(ps1.getMax(1),  12);
    BOOST_CHECK_EQUAL(ps1.getMin(1), -12);
    BOOST_CHECK_CLOSE(ps1.getIntegral(),1,0.1);

    // copy constructor
    vfps::PhaseSpace ps2(ps1);
    BOOST_CHECK_EQUAL(ps2.getMax(0),  12);
    BOOST_CHECK_EQUAL(ps2.getMin(0), -12);
    BOOST_CHECK_EQUAL(ps2.getMax(1),  12);
    BOOST_CHECK_EQUAL(ps2.getMin(1), -12);
    BOOST_CHECK_CLOSE(ps2.getIntegral(),1,0.1f);
    }
    {
    auto nb = 2U;
    // two buckets with 32 grid points for each axis
    vfps::PhaseSpace::resetSize(32,nb);

    vfps::PhaseSpace ps1(-12,12,2,-12,12,4,nullptr,1,1,{{0.5,0.5}});
    BOOST_CHECK_CLOSE(ps1.getIntegral(),1,0.1f);

    for (auto i=0U; i<nb; i++) {
        BOOST_CHECK_CLOSE(ps1.getBunchPopulation()[i],
                          ps1.getSetBunchPopulation()[i],
                          static_cast<vfps::integral_t>(0.1f));
    }
    }
    {
    auto nb = 5U;
    // three buckets with 32 grid points for each axis
    vfps::PhaseSpace::resetSize(32,nb);

    vfps::PhaseSpace ps1(-12,12,2,-12,12,4,nullptr,1,1,{{0.75,0,0.15,0.1,0}});
    BOOST_CHECK_CLOSE(ps1.getIntegral(),1,0.1f);

    for (auto i=0U; i<nb; i++) {
        BOOST_CHECK_CLOSE(ps1.getBunchPopulation()[i],
                          ps1.getSetBunchPopulation()[i],
                          static_cast<vfps::integral_t>(0.1f));
    }
    }
}
