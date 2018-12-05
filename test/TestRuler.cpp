#include <boost/test/unit_test.hpp>

#include <defines.hpp>
#include <PS/Ruler.hpp>

BOOST_AUTO_TEST_CASE( ruler_constructors ){
    vfps::Ruler<vfps::meshaxis_t> axis(15,-4,4);
    BOOST_TEST(axis.zerobin() == 7);
}
