#include <boost/test/unit_test.hpp>

#include <array>

#include "defines.hpp"
#define INOVESA_ALLOW_PS_RESET 1
#include "PS/PhaseSpace.hpp"

BOOST_AUTO_TEST_SUITE( PhaseSpace )

BOOST_AUTO_TEST_CASE( one_bucket ){
    vfps::PhaseSpace::resetSize(32,1);
    // initializing constructor (default data)
    vfps::PhaseSpace ps1(-12,12,3,-12,12,4,nullptr,1,1);
    BOOST_CHECK_EQUAL(ps1.getMax(0),  12);
    BOOST_CHECK_EQUAL(ps1.getMin(0), -12);
    BOOST_CHECK_EQUAL(ps1.getScale(0, "Meter"), 3);
    BOOST_CHECK_EQUAL(ps1.getMax(1),  12);
    BOOST_CHECK_EQUAL(ps1.getMin(1), -12);
    BOOST_CHECK_EQUAL(ps1.getScale(1, "ElectronVolt"), 4);
    BOOST_CHECK_CLOSE(ps1.getIntegral(),1,0.1);

    // copy constructor
    vfps::PhaseSpace ps2(ps1);
    BOOST_CHECK_EQUAL(ps2.getMax(0),  12);
    BOOST_CHECK_EQUAL(ps2.getMin(0), -12);
    BOOST_CHECK_EQUAL(ps2.getMax(1),  12);
    BOOST_CHECK_EQUAL(ps2.getMin(1), -12);
    BOOST_CHECK_CLOSE(ps2.getIntegral(),1,0.1f);


    BOOST_CHECK_CLOSE(ps2.getDelta(0), 24/(32-1.0f), 0.1f);
    BOOST_CHECK_CLOSE(ps2.getDelta(1), 24/(32-1.0f), 0.1f);
}

BOOST_AUTO_TEST_CASE( two_buckets ){
    std::vector<vfps::integral_t> buckets{{0.5,0.5}};
    vfps::PhaseSpace::resetSize(32,buckets.size());

    vfps::PhaseSpace ps1(-12,12,2,-12,12,4,nullptr,1,1,buckets);
    BOOST_CHECK_CLOSE(ps1.getIntegral(),1,0.1f);

    for (auto i=0U; i<buckets.size(); i++) {
        BOOST_CHECK_CLOSE(ps1.getBunchPopulation()[i],
                          ps1.getSetBunchPopulation()[i],
                          static_cast<vfps::integral_t>(0.1f));
    }
}

BOOST_AUTO_TEST_CASE( five_buckets ){
    std::vector<vfps::integral_t> buckets{{0.75,0,0.15,0.1,0}};
    vfps::PhaseSpace::resetSize(32,buckets.size());

    vfps::PhaseSpace ps1(-12,12,2,-12,12,4,nullptr,1,1,buckets);
    BOOST_CHECK_CLOSE(ps1.getIntegral(),1,0.1f);

    for (auto i=0U; i<buckets.size(); i++) {
        BOOST_CHECK_CLOSE(ps1.getBunchPopulation()[i],buckets[i],0.1);
        BOOST_CHECK_CLOSE(ps1.getSetBunchPopulation()[i],buckets[i],0.1);
    }

    ps1.variance(0);
    ps1.variance(1);
    auto mq = ps1.getMoment(0,0);
    auto sq = ps1.getBunchLength();
    auto mp = ps1.getMoment(1,0);
    auto sp = ps1.getEnergySpread();

    for (auto i=0U; i<buckets.size(); i++) {
        BOOST_CHECK_SMALL(mq[i],static_cast<vfps::meshaxis_t>(1e-3));
        BOOST_CHECK_SMALL(mp[i],static_cast<vfps::meshaxis_t>(1e-3));
        if (buckets[i] > 0) {
            BOOST_CHECK_CLOSE(sq[i],1,0.1);
            BOOST_CHECK_CLOSE(sp[i],1,0.1);
        } else {
            BOOST_CHECK_EQUAL(sq[i],0);
            BOOST_CHECK_EQUAL(sp[i],0);
        }
    }
}

BOOST_AUTO_TEST_CASE( no_norm ){
    std::vector<vfps::integral_t> buckets{{0.75,0.75}};
    vfps::PhaseSpace::resetSize(32,buckets.size());

    BOOST_CHECK_THROW(
            vfps::PhaseSpace(-12,12,2,-12,12,4,nullptr,1,1,buckets),
            std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( datacopy ){
    std::vector<vfps::integral_t> buckets{{0.6,0.4}};
    vfps::PhaseSpace::resetSize(2,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.7, 0.7,
                                          0.5, 0.5,

                                          0.5, 0.5,
                                          0.3, 0.3}};

    vfps::PhaseSpace ps1(-1,1,2,-1,1,4,nullptr,1,1,buckets,1,data1.data());

    BOOST_CHECK_EQUAL_COLLECTIONS(data1.data(),data1.data()+data1.size(),
                                  ps1.getData(),ps1.getData()+data1.size());
}

BOOST_AUTO_TEST_CASE( swap ){
    std::vector<vfps::integral_t> buckets{{0.6,0.4}};
    vfps::PhaseSpace::resetSize(2,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.7, 0.7,
                                          0.5, 0.5,

                                          0.5, 0.5,
                                          0.3, 0.3}};

    vfps::PhaseSpace ps1(-1,1,2,-1,1,4,nullptr,1,1,buckets,1,data1.data());


    std::vector<vfps::meshdata_t> data2{{ 0.6, 0.0,
                                          0.0, 0.6,

                                          0.0, 0.4,
                                          0.4, 0.0}};

    vfps::PhaseSpace ps2(-1,1,2,-1,1,4,nullptr,1,1,buckets,1,data2.data());

    BOOST_CHECK_EQUAL_COLLECTIONS(data1.data(),data1.data()+data1.size(),
                                  ps1.getData(),ps1.getData()+ps1.nxyb);
    BOOST_CHECK_EQUAL_COLLECTIONS(data2.data(),data2.data()+data2.size(),
                                  ps2.getData(),ps2.getData()+ps2.nxyb);

    // explicitly use std::swap to check template specialization
    std::swap(ps1,ps2);
    BOOST_CHECK_EQUAL_COLLECTIONS(data1.data(),data1.data()+data1.size(),
                                  ps2.getData(),ps2.getData()+data1.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(data2.data(),data2.data()+data2.size(),
                                  ps1.getData(),ps1.getData()+data2.size());
}

BOOST_AUTO_TEST_CASE( assign ){
    std::vector<vfps::integral_t> buckets{{0.6,0.4}};
    vfps::PhaseSpace::resetSize(2,buckets.size());

    std::vector<vfps::meshdata_t> data1{{ 0.7, 0.7,
                                          0.5, 0.5,

                                          0.5, 0.5,
                                          0.3, 0.3}};

    vfps::PhaseSpace ps1(-1,1,2,-1,1,4,nullptr,1,1,buckets,1,data1.data());

    std::vector<vfps::meshdata_t> data2{{ 0.6, 0.0,
                                          0.0, 0.6,

                                          0.0, 0.4,
                                          0.4, 0.0}};

    vfps::PhaseSpace ps2(-1,1,2,-1,1,4,nullptr,1,1,buckets,1,data2.data());

    BOOST_CHECK_EQUAL_COLLECTIONS(data1.data(),data1.data()+data1.size(),
                                  ps1.getData(),ps1.getData()+ps1.nxyb);
    BOOST_CHECK_EQUAL_COLLECTIONS(data2.data(),data2.data()+data2.size(),
                                  ps2.getData(),ps2.getData()+ps2.nxyb);

    ps2 = ps1;
    BOOST_CHECK_EQUAL_COLLECTIONS(data1.data(),data1.data()+data1.size(),
                                  ps1.getData(),ps1.getData()+ps1.nxyb);
    BOOST_CHECK_EQUAL_COLLECTIONS(data1.data(),data1.data()+data1.size(),
                                  ps2.getData(),ps2.getData()+ps2.nxyb);
}

BOOST_AUTO_TEST_SUITE_END()
