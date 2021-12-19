#include "test_curves.h"

CPPUNIT_TEST_SUITE_REGISTRATION(CurvesTestCase);

void CurvesTestCase::setUp()
{
    fixture1.first.first= "Curve one";
    fixture1.first.second= 1;
    fixture1.second.push_back(1.0);
    fixture1.second.push_back(1.0);

    fixture2.first.first= "Curve two";
    fixture2.first.second= 2;
    fixture2.second.push_back(4.0);
    fixture2.second.push_back(5.0);
}

void CurvesTestCase::tearDown()
{
    fixture1.first.first.clear();
    fixture1.first.second= 0;
    fixture1.second.clear();

    fixture2.first.first.clear();
    fixture2.first.second= 0;
    fixture2.second.clear();
}

void CurvesTestCase::distanceTest()
{
    double distance;
    distance= calculate_distance(fixture1.second, fixture2.second);
    CPPUNIT_ASSERT(distance == 5);
}

void CurvesTestCase::dotProductTest()
{
    int dot_product;
    vector<int> d{4,5};
    dot_product= calculate_dot_product(fixture1, d);
    CPPUNIT_ASSERT( dot_product== 9);
}

void CurvesTestCase::addVectorsTest()
{
    vector <double> sum, result{5.0, 6.0};
    sum = add_vectors(fixture1.second, fixture2.second);
    CPPUNIT_ASSERT(sum == result);
}

void CurvesTestCase::nonZeroCoordinatesTest()
{
    bool non_zero= non_zero_coordinates(fixture1.second);
    CPPUNIT_ASSERT(non_zero == true);
}

void CurvesTestCase::zeroCoordinatesTest()
{
    vector<double> zero_vector{0.0, 0.0};
    bool non_zero= non_zero_coordinates(zero_vector);
    CPPUNIT_ASSERT(non_zero == false);
}

void CurvesTestCase::meanVectorTest()
{
    int i=10;
    vector <double> v{9.0, 6.0, 3.0}, result{3.0, 2.0, 1.0};
    pair<pair<string, int>, vector<double>> mean;
   
    mean=  get_mean_curve_vector(v, 3, i);
    CPPUNIT_ASSERT(mean.second == result);

}