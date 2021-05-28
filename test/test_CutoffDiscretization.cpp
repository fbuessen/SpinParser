#define BOOST_TEST_MODULE "CutoffDiscretizationTest"
#include <boost/test/included/unit_test.hpp>
#include "CutoffDiscretization.hpp"

BOOST_TEST_DONT_PRINT_LOG_VALUE(CutoffIterator);

struct CutoffDiscretizationFixture
{
	CutoffDiscretizationFixture()
	{
		std::vector<float> values({ 5.0, 4.0, 3.0, 2.0, 1.0 });
		c = new CutoffDiscretization(values);
	}

	~CutoffDiscretizationFixture()
	{
		delete c;
	}

	CutoffDiscretization *c;
};

BOOST_FIXTURE_TEST_SUITE(CutoffDiscretizationTest, CutoffDiscretizationFixture);

BOOST_AUTO_TEST_CASE(iterate)
{
	CutoffIterator i = c->begin();
	BOOST_CHECK_EQUAL(*i, 5.0f);

	++i;
	BOOST_CHECK_EQUAL(*i, 4.0f);

	i = c->last();
	BOOST_CHECK_EQUAL(*i, 1.0f);

	++i;
	BOOST_CHECK_EQUAL(i, c->end());
}


BOOST_AUTO_TEST_CASE(find)
{
	CutoffIterator i = c->find(4.0f);
	BOOST_CHECK_EQUAL(*i, 4.0f);

	i = c->find(0.1f);
	BOOST_CHECK_EQUAL(i, c->end());
}

BOOST_AUTO_TEST_SUITE_END();