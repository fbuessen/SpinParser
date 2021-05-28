#define BOOST_TEST_MODULE "IntegratorTest"
#include <boost/test/included/unit_test.hpp>
#include "lib/Integrator.hpp"
#include "lib/ValueBundle.hpp"

class SpinParser
{
public:
	SpinParser(FrequencyDiscretization *f)
	{
		FrgCommon::_frequency = f;
	}

	~SpinParser()
	{
		delete FrgCommon::_frequency;
	}
};

struct IntegratorFixture
{
	IntegratorFixture()
	{
		std::vector<float> values({ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 });
		FrequencyDiscretization *f = new FrequencyDiscretization(values);
		spinParser = new SpinParser(f);
	}

	~IntegratorFixture()
	{
		delete spinParser;
	}

	SpinParser *spinParser;
};

BOOST_FIXTURE_TEST_SUITE(IntegratorTest, IntegratorFixture);

BOOST_AUTO_TEST_CASE(IntegratorIntegrate, *boost::unit_test::tolerance(1.0e-6))
{
	std::function<float(float)> integrandLinear = [](float x)->float { return x; };
	std::function<float(float)> integrandQuadratic = [](float x)->float { return x*x; };
	auto min = FrgCommon::frequency().begin();
	auto max = FrgCommon::frequency().last();

	BOOST_CHECK_EQUAL(Integrator::integrate(min, max, integrandLinear), 49.5f);
	BOOST_CHECK_EQUAL(Integrator::integrate(min, max, integrandQuadratic), 334.5f);
}

BOOST_AUTO_TEST_CASE(IntegratorIntegrateWithObscureLeftBoundary, *boost::unit_test::tolerance(1.0e-6))
{
	std::function<float(float)> integrandLinear = [](float x)->float { return x; };
	std::function<float(float)> integrandQuadratic = [](float x)->float { return x * x; };
	auto min = 1.5f;
	auto max = FrgCommon::frequency().last();

	BOOST_CHECK_EQUAL(Integrator::integrateWithObscureLeftBoundary(min, max, integrandLinear), 48.875f);
	BOOST_CHECK_EQUAL(Integrator::integrateWithObscureLeftBoundary(min, max, integrandQuadratic), 333.5625f);
}

BOOST_AUTO_TEST_CASE(IntegratorIntegrateWithObscureRightBoundary, *boost::unit_test::tolerance(1.0e-6))
{
	std::function<float(float)> integrandLinear = [](float x)->float { return x; };
	std::function<float(float)> integrandQuadratic = [](float x)->float { return x * x; };
	auto min = FrgCommon::frequency().begin();
	auto max = 9.5f;

	BOOST_CHECK_EQUAL(Integrator::integrateWithObscureRightBoundary(min, max, integrandLinear), 44.625f);
	BOOST_CHECK_EQUAL(Integrator::integrateWithObscureRightBoundary(min, max, integrandQuadratic), 286.8125f);
}

BOOST_AUTO_TEST_CASE(IntegratorIntegrateWithObscureBoundaries, *boost::unit_test::tolerance(1.0e-6))
{
	std::function<float(float)> integrandLinear = [](float x)->float { return x; };
	std::function<float(float)> integrandQuadratic = [](float x)->float { return x * x; };
	auto min = 1.5f;
	auto max = 9.5f;

	BOOST_CHECK_EQUAL(Integrator::integrateWithObscureBoundaries(min, max, integrandLinear), 44.0f);
	BOOST_CHECK_EQUAL(Integrator::integrateWithObscureBoundaries(min, max, integrandQuadratic), 285.875f);
}

BOOST_AUTO_TEST_CASE(ImplicitIntegratorIntegrateWithObscureLeftBoundary, *boost::unit_test::tolerance(1.0e-6))
{
	std::function<void(float, ValueSuperbundle<float, 1> &)> integrand = [](float x, ValueSuperbundle<float, 1> &out)->void { out.bundle(0)[0] = x; out.bundle(0)[1] = x * x; };
	auto min = 1.5f;
	auto max = FrgCommon::frequency().last();

	ValueSuperbundle<float, 1> buffer = ValueSuperbundle<float, 1>(2);
	ValueSuperbundle<float, 1> result = ValueSuperbundle<float, 1>(2);
	ImplicitIntegrator::integrateWithObscureLeftBoundary<ValueSuperbundle<float, 1> >(min, max, integrand, buffer, result);
	BOOST_CHECK_EQUAL(result.bundle(0)[0], 48.875f);
	BOOST_CHECK_EQUAL(result.bundle(0)[1], 333.5625f);
}

BOOST_AUTO_TEST_CASE(ImplicitIntegratorIntegrateWithObscureRightBoundary, *boost::unit_test::tolerance(1.0e-6))
{
	std::function<void(float, ValueSuperbundle<float, 1> &)> integrand = [](float x, ValueSuperbundle<float, 1> &out)->void { out.bundle(0)[0] = x; out.bundle(0)[1] = x * x; };
	auto min = FrgCommon::frequency().begin();
	auto max = 9.5f;

	ValueSuperbundle<float, 1> buffer = ValueSuperbundle<float, 1>(2);
	ValueSuperbundle<float, 1> result = ValueSuperbundle<float, 1>(2);
	ImplicitIntegrator::integrateWithObscureRightBoundary<ValueSuperbundle<float, 1> >(min, max, integrand, buffer, result);
	BOOST_CHECK_EQUAL(result.bundle(0)[0], 44.625f);
	BOOST_CHECK_EQUAL(result.bundle(0)[1], 286.8125f);
}

BOOST_AUTO_TEST_CASE(ImplicitIntegratorIntegrateWithObscureBoundaries, *boost::unit_test::tolerance(1.0e-6))
{
	std::function<void(float, ValueSuperbundle<float, 1> &)> integrand = [](float x, ValueSuperbundle<float, 1> &out)->void { out.bundle(0)[0] = x; out.bundle(0)[1] = x * x; };
	auto min = 1.5f;
	auto max = 9.5f;

	ValueSuperbundle<float, 1> buffer = ValueSuperbundle<float, 1>(2);
	ValueSuperbundle<float, 1> result = ValueSuperbundle<float, 1>(2);
	ImplicitIntegrator::integrateWithObscureBoundaries<ValueSuperbundle<float, 1> >(min, max, integrand, buffer, result);
	BOOST_CHECK_EQUAL(result.bundle(0)[0], 44.0f);
	BOOST_CHECK_EQUAL(result.bundle(0)[1], 285.875f);
}

BOOST_AUTO_TEST_SUITE_END();