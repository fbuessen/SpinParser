#define BOOST_TEST_MODULE "XYZVertexTwoParticleTest"
#include <boost/test/included/unit_test.hpp>
#include "LatticeModelFactory.hpp"
#include "XYZ/XYZVertexTwoParticle.hpp"

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

struct XYZVertexTwoParticleFixture
{
	XYZVertexTwoParticleFixture()
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
		v = new XYZVertexTwoParticle;
	}

	~XYZVertexTwoParticleFixture()
	{
		delete spinParser;
		delete v;
	}

	XYZVertexTwoParticle *v;
	SpinParser *spinParser;
};

BOOST_FIXTURE_TEST_SUITE(XYZVertexTwoParticleTest, XYZVertexTwoParticleFixture);

BOOST_AUTO_TEST_CASE(ExpandIteratorFull)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SpinComponent::X) = float(i);
		v->getValueRef(i, SpinComponent::Y) = float(i + 1);
		v->getValueRef(i, SpinComponent::Z) = float(i + 2);
		v->getValueRef(i, SpinComponent::None) = float(i + 3);
	}
	for (int i = 0; i < v->size; ++i)
	{
		LatticeIterator i1;
		float s = 0.0;
		float t = 0.0;
		float u = 0.0;
		v->expandIterator(i, i1, s, t, u);

		BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), i1, s, t, u, SpinComponent::X, XYZVertexTwoParticle::FrequencyChannel::All), float(i));
		BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), i1, s, t, u, SpinComponent::Y, XYZVertexTwoParticle::FrequencyChannel::All), float(i + 1));
		BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), i1, s, t, u, SpinComponent::Z, XYZVertexTwoParticle::FrequencyChannel::All), float(i + 2));
		BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), i1, s, t, u, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All), float(i + 3));
	}
}

BOOST_AUTO_TEST_CASE(ExpandIteratorFrequency)
{
	for (int i = 0; i < v->sizeFrequency; ++i)
	{
		for (int rid = 0; rid < FrgCommon::lattice().size; ++rid)
		{
			v->getValueRef(i * FrgCommon::lattice().size + rid, SpinComponent::X) = float(i);
			v->getValueRef(i * FrgCommon::lattice().size + rid, SpinComponent::Y) = float(i + 1);
			v->getValueRef(i * FrgCommon::lattice().size + rid, SpinComponent::Z) = float(i + 2);
			v->getValueRef(i * FrgCommon::lattice().size + rid, SpinComponent::None) = float(i + 3);
		}
	}
	for (int i = 0; i < v->sizeFrequency; ++i)
	{
		for (int rid = 0; rid < FrgCommon::lattice().size; ++rid)
		{
			float s = 0.0;
			float t = 0.0;
			float u = 0.0;
			v->expandIterator(i, s, t, u);

			BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().fromParametrization(rid), s, t, u, SpinComponent::X, XYZVertexTwoParticle::FrequencyChannel::All), float(i));
			BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().fromParametrization(rid), s, t, u, SpinComponent::Y, XYZVertexTwoParticle::FrequencyChannel::All), float(i + 1));
			BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().fromParametrization(rid), s, t, u, SpinComponent::Z, XYZVertexTwoParticle::FrequencyChannel::All), float(i + 2));
			BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().fromParametrization(rid), s, t, u, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All), float(i + 3));
		}
	}
}

BOOST_AUTO_TEST_CASE(getValueSymmetry)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SpinComponent::X) = float(i);
		v->getValueRef(i, SpinComponent::Y) = float(i + 1);
		v->getValueRef(i, SpinComponent::Z) = float(i + 2);
		v->getValueRef(i, SpinComponent::None) = float(i + 3);
	}

	SublatticeIterator i1 = FrgCommon::lattice().getBasis();
	SublatticeIterator i2 = ++FrgCommon::lattice().getRange(i1);
	FrequencyIterator s = FrgCommon::frequency().begin();
	FrequencyIterator t = ++FrgCommon::frequency().begin();
	FrequencyIterator u = ++++FrgCommon::frequency().begin();

	for (int i = 0; i < 3; ++i)
	{
		BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All), v->getValue(i2, i1, -*s, *t, *u, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All));
		BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All), v->getValue(i1, i2, *s, -*t, *u, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All));
		BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All), v->getValue(i2, i1, *s, *t, -*u, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All));
		BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All), v->getValue(i1, i2, *u, *t, *s, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All));
	}

	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All), v->getValue(i2, i1, -*s, *t, *u, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All));
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All), v->getValue(i1, i2, *s, -*t, *u, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All));
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All), v->getValue(i2, i1, *s, *t, -*u, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All));
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All), -v->getValue(i1, i2, *u, *t, *s, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::All));
}

BOOST_AUTO_TEST_CASE(getValueInterpolation)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SpinComponent::X) = float(i);
		v->getValueRef(i, SpinComponent::Y) = float(i + 1);
		v->getValueRef(i, SpinComponent::Z) = float(i + 2);
		v->getValueRef(i, SpinComponent::None) = float(i + 3);
	}

	SublatticeIterator i1 = FrgCommon::lattice().getBasis();
	SublatticeIterator i2 = ++FrgCommon::lattice().getRange(i1);
	FrequencyIterator wExact = FrgCommon::frequency().begin();
	FrequencyIterator wLower = ++FrgCommon::frequency().begin();
	FrequencyIterator wUpper= ++++FrgCommon::frequency().begin();
	float lambda = 0.2f;
	float vLower, vUpper, vExact;
	
	for (int i = 0; i < 4; ++i)
	{
		vLower = v->getValue(i1, i2, *wLower, *wExact, *wExact, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All);
		vUpper = v->getValue(i1, i2, *wUpper, *wExact, *wExact, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All);
		vExact = v->getValue(i1, i2, (1.0f - lambda) * *wLower + lambda * *wUpper, *wExact, *wExact, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::None);
		BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);

		vLower = v->getValue(i1, i2, *wExact, *wLower, *wExact, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All);
		vUpper = v->getValue(i1, i2, *wExact, *wUpper, *wExact, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All);
		vExact = v->getValue(i1, i2, *wExact, (1.0f - lambda) * *wLower + lambda * *wUpper, *wExact, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::S);
		BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);

		vLower = v->getValue(i1, i2, *wExact, *wExact, *wLower, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All);
		vUpper = v->getValue(i1, i2, *wExact, *wExact, *wUpper, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All);
		vExact = v->getValue(i1, i2, *wExact, *wExact, (1.0f - lambda) * *wLower + lambda * *wUpper, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::T);
		BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);

		vLower = v->getValue(i1, i2, *wLower, *wExact, *wExact, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All);
		vUpper = v->getValue(i1, i2, *wUpper, *wExact, *wExact, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All);
		vExact = v->getValue(i1, i2, (1.0f - lambda) * *wLower + lambda * *wUpper, *wExact, *wExact, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::U);
		BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);
	}
}

BOOST_AUTO_TEST_CASE(getValueBuffered)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SpinComponent::X) = float(i);
		v->getValueRef(i, SpinComponent::Y) = float(i + 1);
		v->getValueRef(i, SpinComponent::Z) = float(i + 2);
		v->getValueRef(i, SpinComponent::None) = float(i + 3);
	}

	SublatticeIterator i1 = FrgCommon::lattice().getBasis();
	SublatticeIterator i2 = ++FrgCommon::lattice().getRange(i1);

	for (int i = 0; i < 4; ++i)
	{
		auto ab = v->generateAccessBuffer(1.1f, 2.2f, 3.3f);
		BOOST_CHECK_CLOSE(v->getValue(i1, i2, static_cast<SpinComponent>(i), ab), v->getValue(i1, i2, 1.1f, 2.2f, 3.3f, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::None), 0.0001);

		ab = v->generateAccessBuffer(2.0f, 3.0f, 4.0f);
		BOOST_CHECK_CLOSE(v->getValue(i1, i2, static_cast<SpinComponent>(i), ab), v->getValue(i1, i2, 2.0f, 3.0f, 4.0f, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::All), 0.0001);

		auto ab2 = v->generateAccessBuffer(2.0f, 3.5f, 4.5f, XYZVertexTwoParticle::FrequencyChannel::S);
		BOOST_CHECK_CLOSE(v->getValue(i1, i2, static_cast<SpinComponent>(i), ab2), v->getValue(i1, i2, 2.0f, 3.5f, 4.5f, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::S), 0.0001);

		ab2 = v->generateAccessBuffer(2.5f, 3.0f, 4.5f, XYZVertexTwoParticle::FrequencyChannel::T);
		BOOST_CHECK_CLOSE(v->getValue(i1, i2, static_cast<SpinComponent>(i), ab2), v->getValue(i1, i2, 2.5f, 3.0f, 4.5f, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::T), 0.0001);

		ab2 = v->generateAccessBuffer(2.5f, 3.5f, 4.0f, XYZVertexTwoParticle::FrequencyChannel::U);
		BOOST_CHECK_CLOSE(v->getValue(i1, i2, static_cast<SpinComponent>(i), ab2), v->getValue(i1, i2, 2.5f, 3.5f, 4.0f, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::U), 0.0001);
	}
}

BOOST_AUTO_TEST_CASE(getValueLocal)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SpinComponent::X) = float(i);
		v->getValueRef(i, SpinComponent::Y) = float(i + 1);
		v->getValueRef(i, SpinComponent::Z) = float(i + 2);
		v->getValueRef(i, SpinComponent::None) = float(i + 3);
	}

	for (int i = 0; i < 4; ++i)
	{
		auto ab = v->generateAccessBuffer(1.1f, 2.2f, 3.3f);
		BOOST_CHECK_CLOSE(v->getValueLocal(static_cast<SpinComponent>(i), ab), v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), 1.1f, 2.2f, 3.3f, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::None), 0.0001);
	}
}

BOOST_AUTO_TEST_CASE(getValueSuperbundle)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SpinComponent::X) = float(i);
		v->getValueRef(i, SpinComponent::Y) = float(i + 1);
		v->getValueRef(i, SpinComponent::Z) = float(i + 2);
		v->getValueRef(i, SpinComponent::None) = float(i + 3);
	}

	auto ab = v->generateAccessBuffer(1.1f, 2.2f, 3.3f);
	ValueSuperbundle<float, 4> b(FrgCommon::lattice().size);
	v->getValueSuperbundle(ab, b);

	for (int rid = 0; rid < FrgCommon::lattice().size; ++rid)
	{
		for (int i = 0; i < 4; ++i)
		{
			BOOST_CHECK_CLOSE(b.bundle(i)[rid], v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().fromParametrization(rid), 1.1f, 2.2f, 3.3f, static_cast<SpinComponent>(i), XYZVertexTwoParticle::FrequencyChannel::None), 0.0001);
		}
	}
}

BOOST_AUTO_TEST_SUITE_END();