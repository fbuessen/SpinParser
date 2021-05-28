#define BOOST_TEST_MODULE "SU2VertexSingleParticleTest"
#include <boost/test/included/unit_test.hpp>
#include "LatticeModelFactory.hpp"
#include "SU2/SU2VertexSingleParticle.hpp"

class SpinParser
{
public:
	SpinParser(Lattice *l, FrequencyDiscretization *f)
	{
		FrgCommon::_lattice = l;
		FrgCommon::_frequency = f;
	}

	~SpinParser()
	{
		delete FrgCommon::_lattice;
		delete FrgCommon::_frequency;
	}
};

struct SU2VertexSingleParticleFixture
{
	SU2VertexSingleParticleFixture()
	{
		//construct square lattice
		LatticeModelFactory::LatticeUnitCell uc;
		uc.basisSites.push_back(geometry::Vec3<double>(0.0, 0.0, 0.0));
		uc.latticeVectors.push_back(geometry::Vec3<double>(1.0, 0.0, 0.0));
		uc.latticeVectors.push_back(geometry::Vec3<double>(0.0, 1.0, 0.0));
		uc.latticeVectors.push_back(geometry::Vec3<double>(0.0, 0.0, 1.0));
		uc.latticeBonds.push_back(LatticeModelFactory::LatticeBond(0, 0, 1, 0, 0));
		uc.latticeBonds.push_back(LatticeModelFactory::LatticeBond(0, 0, 0, 1, 0));

		LatticeModelFactory::SpinModelUnitCell model;
		LatticeModelFactory::SpinInteraction i1(LatticeModelFactory::LatticeSite(0, 0, 0, 0), LatticeModelFactory::LatticeSite(1, 0, 0, 0));
		i1.interactionStrength[0][0] = 1.0f;
		i1.interactionStrength[1][1] = 1.0f;
		i1.interactionStrength[2][2] = 1.0f;
		model.interactions.push_back(i1);
		LatticeModelFactory::SpinInteraction i2(LatticeModelFactory::LatticeSite(0, 0, 0, 0), LatticeModelFactory::LatticeSite(0, 1, 0, 0));
		i2.interactionStrength[0][0] = 1.0f;
		i2.interactionStrength[1][1] = 1.0f;
		i2.interactionStrength[2][2] = 1.0f;
		model.interactions.push_back(i2);

		Log::log << Log::setDisplayLogLevel(Log::LogLevel::None);
		std::pair<Lattice *, SpinModel *> product = LatticeModelFactory::newLatticeModel(uc, model, 3);
		Lattice *l = product.first;
		delete product.second;

		//construct frequency discretization
		std::vector<float> values({ 1.0, 2.0, 3.0, 4.0, 5.0 });
		FrequencyDiscretization *f = new FrequencyDiscretization(values);

		//write lattice and frequency to FrgCommon
		spinParser = new SpinParser(l, f);

		//construct vertex
		v = new SU2VertexSingleParticle;
	}

	~SU2VertexSingleParticleFixture()
	{
		delete spinParser;
		delete v;
	}

	SU2VertexSingleParticle *v;
	SpinParser *spinParser;
};

BOOST_FIXTURE_TEST_SUITE(SU2VertexSingleParticleTest, SU2VertexSingleParticleFixture);

BOOST_AUTO_TEST_CASE(ExpandIterator)
{
	for (int i = 0; i < v->size; ++i) v->getValueRef(i) = float(i);
	for (int i = 0; i < v->size; ++i)
	{
		float w;
		v->expandIterator(i, w);
		BOOST_CHECK_EQUAL(v->getValue(w), float(i));
	}
}

BOOST_AUTO_TEST_CASE(getValueSymmetry)
{
	for (int i = 0; i < v->size; ++i) v->getValueRef(i) = float(i);
	for (auto w = FrgCommon::frequency().begin(); w != FrgCommon::frequency().end(); ++w)
	{
		BOOST_CHECK_EQUAL(v->getValue(*w), -v->getValue(-*w));
	}
}

BOOST_AUTO_TEST_CASE(getValueInterpolation)
{
	float lower = 2.0f;
	float upper = 3.0f;
	v->getValueRef(0) = lower;
	v->getValueRef(1) = upper;

	float lambda = 0.2f;
	float w = (1.0f - lambda) * *FrgCommon::frequency().begin() + lambda * *(++FrgCommon::frequency().begin());
	BOOST_CHECK_EQUAL(v->getValue(w), (1.0f - lambda) * lower + lambda * upper);
}

BOOST_AUTO_TEST_SUITE_END();