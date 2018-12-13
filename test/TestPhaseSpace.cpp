#include <boost/test/unit_test.hpp>

#include <defines.hpp>
#include <PS/PhaseSpace.hpp>

BOOST_AUTO_TEST_CASE( phasespace_constructors ){
    vfps::PhaseSpace::setSize(32,32,1);
    vfps::PhaseSpace ps(32,-1,1,3,-2,2,4,nullptr, 1,1);
    BOOST_CHECK(ps.getMax(0) ==  1);
    BOOST_CHECK(ps.getMin(0) == -1);
    BOOST_CHECK(ps.getMax(1) ==  2);
    BOOST_CHECK(ps.getMin(1) == -2);
}
