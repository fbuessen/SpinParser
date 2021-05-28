#define BOOST_TEST_MODULE "FrequencyDiscretizationTest"
#include <boost/test/included/unit_test.hpp>
#include "FrequencyDiscretization.hpp"

BOOST_TEST_DONT_PRINT_LOG_VALUE(FrequencyIterator);

struct FrequencyDiscretizationFixture
{
	FrequencyDiscretizationFixture()
	{
		std::vector<float> values({ 1.0, 2.0, 3.0, 4.0, 5.0 });
		f = new FrequencyDiscretization(values);
	}

	~FrequencyDiscretizationFixture()
	{
		delete f;
	}

	FrequencyDiscretization* f;
};

BOOST_FIXTURE_TEST_SUITE(FrequencyDiscretizationTest, FrequencyDiscretizationFixture);

BOOST_AUTO_TEST_CASE(iterate)
{
	FrequencyIterator i = f->begin();
	BOOST_CHECK_EQUAL(*i, 1.0f);

	++i;
	BOOST_CHECK_EQUAL(*i, 2.0f);

	i = f->last();
	BOOST_CHECK_EQUAL(*i, 5.0f);

	++i;
	BOOST_CHECK_EQUAL(i, f->end());

	i = f->beginNegative();
	BOOST_CHECK_EQUAL(*i, -5.0f);

	++i;
	BOOST_CHECK_EQUAL(*i, -4.0f);
}

BOOST_AUTO_TEST_CASE(lesser)
{
	FrequencyIterator i = f->lesser(10.0f);
	BOOST_CHECK_EQUAL(*i, 5.0f);

	i = f->lesser(0.5f);
	BOOST_CHECK_EQUAL(*i, 1.0f);
	
	i = f->lesser(3.5f);
	BOOST_CHECK_EQUAL(*i, 3.0f);

	i = f->lesser(-10.0f);
	BOOST_CHECK_EQUAL(*i, -5.0f);

	i = f->lesser(-0.5f);
	BOOST_CHECK_EQUAL(*i, -1.0f);

	i = f->lesser(-3.5f);
	BOOST_CHECK_EQUAL(*i, -4.0f);
}

BOOST_AUTO_TEST_CASE(greater)
{
	FrequencyIterator i = f->greater(10.0f);
	BOOST_CHECK_EQUAL(*i, 5.0f);

	i = f->greater(0.5f);
	BOOST_CHECK_EQUAL(*i, 1.0f);

	i = f->greater(3.5f);
	BOOST_CHECK_EQUAL(*i, 4.0f);

	i = f->greater(-10.0f);
	BOOST_CHECK_EQUAL(*i, -5.0f);

	i = f->greater(-0.5f);
	BOOST_CHECK_EQUAL(*i, -1.0f);

	i = f->greater(-3.5f);
	BOOST_CHECK_EQUAL(*i, -3.0f);
}

BOOST_AUTO_TEST_CASE(offset)
{
	int offset = f->offset(1.0f);
	BOOST_CHECK_EQUAL(offset, 0);

	offset = f->offset(2.0f);
	BOOST_CHECK_EQUAL(offset, 1);

	offset = f->offset(5.0f);
	BOOST_CHECK_EQUAL(offset, 4);
}

BOOST_AUTO_TEST_CASE(interpolateOffset)
{
	int lowerOffset, upperOffset;
	float bias;
	f->interpolateOffset(10.0f, lowerOffset, upperOffset, bias);
	BOOST_CHECK_EQUAL(lowerOffset, 4);
	BOOST_CHECK_EQUAL(upperOffset, 4);
	BOOST_CHECK_GE(bias, 0.0f);
	BOOST_CHECK_LE(bias, 1.0f);

	f->interpolateOffset(4.5f, lowerOffset, upperOffset, bias);
	BOOST_CHECK_EQUAL(lowerOffset, 3);
	BOOST_CHECK_EQUAL(upperOffset, 4);
	BOOST_CHECK_CLOSE(bias, 0.5f, 0.0001);

	f->interpolateOffset(0.1f, lowerOffset, upperOffset, bias);
	BOOST_CHECK_EQUAL(lowerOffset, 0);
	BOOST_CHECK_EQUAL(upperOffset, 0);
	BOOST_CHECK_GE(bias, 0.0f);
	BOOST_CHECK_LE(bias, 1.0f);
}


BOOST_AUTO_TEST_SUITE_END();