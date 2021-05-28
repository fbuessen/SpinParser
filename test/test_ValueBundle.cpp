#define BOOST_TEST_MODULE "ValueBundleTest"
#include <boost/test/included/unit_test.hpp>
#include "lib/ValueBundle.hpp"


BOOST_AUTO_TEST_SUITE(ValueBundleTest);

BOOST_AUTO_TEST_CASE(ValueBundleArithmetic)
{
	const int dataSize = 16;
	float data1[dataSize];
	float data2[dataSize];
	float data3[dataSize];

	auto reset = [](ValueBundle<float> &b)->void {
		for (int i = 0; i < dataSize; ++i) b[i] = float(i);
	};

	ValueBundle<float> b1 = ValueBundle<float>(&data1[0], dataSize);
	ValueBundle<float> b2 = ValueBundle<float>(&data2[0], dataSize);
	ValueBundle<float> b3 = ValueBundle<float>(&data3[0], dataSize);

	BOOST_CHECK_EQUAL(b1.size(), dataSize);
	
	reset(b1);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i));

	reset(b1);
	reset(b2);
	b1.multAdd(2.0f, b2);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i + 2.0f * i));

	reset(b1);
	reset(b2);
	b1.multAdd(b2, 2.0f);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i + 2.0f * i));

	reset(b1);
	reset(b2);
	reset(b3);
	b1.multAdd(b2, b3);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i + i * i));

	reset(b1);
	reset(b2);
	reset(b3);
	b1.multAdd(2.0f, b2, b3);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i + 2.0f * i * i));

	reset(b1);
	reset(b2);
	b1.multSub(2.0f, b2);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i - 2.0f * i));

	reset(b1);
	reset(b2);
	b1.multSub(b2, 2.0f);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i - 2.0f * i));

	reset(b1);
	reset(b2);
	reset(b3);
	b1.multSub(b2, b3);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i - i * i));

	reset(b1);
	reset(b2);
	reset(b3);
	b1.multSub(2.0f, b2, b3);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i - 2.0f * i * i));

	reset(b1);
	reset(b2);
	b1 += b2;
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i + i));

	reset(b1);
	reset(b2);
	b1 -= b2;
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i - i));

	reset(b1);
	b1 *= 2.0f;
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i * 2.0f));

	reset(b1);
	b1 /= 2.0f;
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(b1[i], float(i / 2.0f));
}

BOOST_AUTO_TEST_CASE(ValueSuperbundleArithmetic)
{
	const int dataSize = 16;

	auto reset = [](ValueSuperbundle<float,2> &s)->void {
		for (int i = 0; i < 2; ++i)
		{
			auto b = s.bundle(i);
			for (int j = 0; j < dataSize; ++j) b[j] = float(j);
		}
	};

	ValueSuperbundle<float, 2> s1 = ValueSuperbundle<float, 2>(dataSize);
	ValueSuperbundle<float, 2> s2 = ValueSuperbundle<float, 2>(dataSize);

	s1.reset();
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(0)[i], 0.0f);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(1)[i], 0.0f);

	ValueSuperbundle<float, 2> s3 = ValueSuperbundle<float, 2>(s1);
	s1.bundle(0)[0] = 42.0f;
	s3.bundle(0)[1] = 43.0f;
	BOOST_CHECK_EQUAL(s3.bundle(0)[0], 42.0f);
	BOOST_CHECK_EQUAL(s1.bundle(0)[1], 43.0f);

	reset(s1);
	reset(s2);
	s1.multAdd(2.0f, s2);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(0)[i], float(i + 2.0f * i));
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(1)[i], float(i + 2.0f * i));

	reset(s1);
	s1 *= 2.0f;
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(0)[i], float(2.0f * i));
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(1)[i], float(2.0f * i));

	reset(s1);
	s1 /= 2.0f;
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(0)[i], float(i) / 2.0f);
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(1)[i], float(i) / 2.0f);

	reset(s1);
	reset(s2);
	s1 += s2;
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(0)[i], float(2.0f * i));
	for (int i = 0; i < dataSize; ++i) BOOST_CHECK_EQUAL(s1.bundle(1)[i], float(2.0f * i));
}
BOOST_AUTO_TEST_SUITE_END();