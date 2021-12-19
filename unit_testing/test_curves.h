#ifndef TEST_CURVES_H
#define TEST_CURVES_H

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <utility>
#include <vector>
#include "vector_ops.h"

class CurvesTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CurvesTestCase);
    CPPUNIT_TEST(distanceTest);
    CPPUNIT_TEST(dotProductTest);
    CPPUNIT_TEST(addVectorsTest);
    CPPUNIT_TEST(nonZeroCoordinatesTest);
    CPPUNIT_TEST(zeroCoordinatesTest);
    CPPUNIT_TEST(meanVectorTest);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void setUp();
    void tearDown();

protected:
    void distanceTest(void);
    void dotProductTest(void);
    void addVectorsTest(void);
    void nonZeroCoordinatesTest(void);
    void zeroCoordinatesTest(void);
    void meanVectorTest(void);

private:
    pair<pair<string, int>, vector<double>> fixture1;
    pair<pair<string, int>, vector<double>> fixture2;
};



#endif