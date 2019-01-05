#include <boost/test/unit_test.hpp>

#include <defines.hpp>
#include <PS/PhaseSpace.hpp>

BOOST_AUTO_TEST_CASE( phasespace_constructors ){
    vfps::PhaseSpace::setSize(32,1);

    // default constructor
    vfps::PhaseSpace ps1(-12,12,3,-22,22,4,nullptr,1,1);
    BOOST_CHECK_EQUAL(ps1.getMax(0),  12);
    BOOST_CHECK_EQUAL(ps1.getMin(0), -12);
    BOOST_CHECK_EQUAL(ps1.getMax(1),  22);
    BOOST_CHECK_EQUAL(ps1.getMin(1), -22);
    BOOST_CHECK_CLOSE(ps1.getIntegral(),
                      static_cast<vfps::integral_t>(1),
                      static_cast<vfps::integral_t>(0.1f));

    // copy constructor
    vfps::PhaseSpace ps2(ps1);
    BOOST_CHECK_EQUAL(ps2.getMax(0),  12);
    BOOST_CHECK_EQUAL(ps2.getMin(0), -12);
    BOOST_CHECK_EQUAL(ps2.getMax(1),  22);
    BOOST_CHECK_EQUAL(ps2.getMin(1), -22);
    BOOST_CHECK_CLOSE(ps2.getIntegral(),
                      static_cast<vfps::integral_t>(1),
                      static_cast<vfps::integral_t>(0.1f));
}
