#include <boost/test/unit_test.hpp>

#include <defines.hpp>
#include <PS/PhaseSpace.hpp>

BOOST_AUTO_TEST_CASE( phasespace_constructors ){
    vfps::PhaseSpace::setSize(32,32,1);
    vfps::PhaseSpace ps(-12,12,3,-22,22,4,nullptr,1,1);
    BOOST_CHECK_EQUAL(ps.getMax(0),  12);
    BOOST_CHECK_EQUAL(ps.getMin(0), -12);
    BOOST_CHECK_EQUAL(ps.getMax(1),  22);
    BOOST_CHECK_EQUAL(ps.getMin(1), -22);
    BOOST_CHECK_CLOSE(ps.getIntegral(),
                      static_cast<vfps::integral_t>(1),
                      static_cast<vfps::integral_t>(0.1f));
}
