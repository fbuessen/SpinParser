#define BOOST_TEST_MODULE "SU2VertexTwoParticleTest"
#include <boost/test/included/unit_test.hpp>
#include "LatticeModelFactory.hpp"
#include "SU2/SU2VertexTwoParticle.hpp"

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

struct SU2VertexTwoParticleFixture
{
	SU2VertexTwoParticleFixture()
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
		v = new SU2VertexTwoParticle;
	}

	~SU2VertexTwoParticleFixture()
	{
		delete spinParser;
		delete v;
	}

	SU2VertexTwoParticle *v;
	SpinParser *spinParser;
};

BOOST_FIXTURE_TEST_SUITE(SU2VertexTwoParticleTest, SU2VertexTwoParticleFixture);

BOOST_AUTO_TEST_CASE(ExpandIteratorFull)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Spin) = float(i);
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Density) = float(i) + 1.0f;
	}
	for (int i = 0; i < v->size; ++i)
	{
		LatticeIterator i1;
		float s = 0.0;
		float t = 0.0;
		float u = 0.0;
		v->expandIterator(i, i1, s, t, u);

		BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), i1, s, t, u, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All), float(i));
		BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), i1, s, t, u, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All), float(i) + 1.0f);
	}
}

BOOST_AUTO_TEST_CASE(ExpandIteratorFrequency)
{
	for (int i = 0; i < v->sizeFrequency; ++i)
	{
		for (int rid = 0; rid < FrgCommon::lattice().size; ++rid)
		{
			v->getValueRef(i * FrgCommon::lattice().size + rid, SU2VertexTwoParticle::Symmetry::Spin) = float(i);
			v->getValueRef(i * FrgCommon::lattice().size + rid, SU2VertexTwoParticle::Symmetry::Density) = float(i) + 1.0f;
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

			BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().fromParametrization(rid), s, t, u, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All), float(i));
			BOOST_CHECK_EQUAL(v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().fromParametrization(rid), s, t, u, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All), float(i) + 1.0f);
		}
	}
}

BOOST_AUTO_TEST_CASE(getValueSymmetry)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Spin) = float(i);
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Density) = float(i) + 1.0f;
	}

	SublatticeIterator i1 = FrgCommon::lattice().getBasis(); 
	SublatticeIterator i2 = ++FrgCommon::lattice().getRange(i1); 
	FrequencyIterator s = FrgCommon::frequency().begin();
	FrequencyIterator t = ++FrgCommon::frequency().begin();
	FrequencyIterator u = ++++FrgCommon::frequency().begin();
		
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All), v->getValue(i2, i1, -*s, *t, *u, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All));
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All), v->getValue(i1, i2, *s, -*t, *u, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All));
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All), v->getValue(i2, i1, *s, *t, -*u, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All));
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All), v->getValue(i1, i2, *u, *t, *s, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All));
						
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All), v->getValue(i2, i1, -*s, *t, *u, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All));
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All), v->getValue(i1, i2, *s, -*t, *u, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All));
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All), v->getValue(i2, i1, *s, *t, -*u, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All));
	BOOST_CHECK_EQUAL(v->getValue(i1, i2, *s, *t, *u, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All), -v->getValue(i1, i2, *u, *t, *s, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All));
}

BOOST_AUTO_TEST_CASE(getValueInterpolation)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Spin) = float(i);
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Density) = float(i) + 1.0f;
	}

	SublatticeIterator i1 = FrgCommon::lattice().getBasis();
	SublatticeIterator i2 = ++FrgCommon::lattice().getRange(i1);
	FrequencyIterator wExact = FrgCommon::frequency().begin();
	FrequencyIterator wLower = ++FrgCommon::frequency().begin();
	FrequencyIterator wUpper= ++++FrgCommon::frequency().begin();
	float lambda = 0.2f;
	float vLower, vUpper, vExact;
	
	vLower = v->getValue(i1, i2,*wLower, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All);
	vUpper= v->getValue(i1, i2, *wUpper, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All);
	vExact = v->getValue(i1, i2, (1.0f - lambda) * *wLower + lambda * *wUpper, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::None);
	BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);

	vLower = v->getValue(i1, i2, *wExact, *wLower, *wExact, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All);
	vUpper = v->getValue(i1, i2, *wExact, *wUpper, *wExact, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All);
	vExact = v->getValue(i1, i2, *wExact, (1.0f - lambda) * *wLower + lambda * *wUpper , *wExact, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::S);
	BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);

	vLower = v->getValue(i1, i2, *wExact, *wExact, *wLower, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All);
	vUpper = v->getValue(i1, i2, *wExact, *wExact, *wUpper, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All);
	vExact = v->getValue(i1, i2, *wExact, *wExact, (1.0f - lambda) * *wLower + lambda * *wUpper, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::T);
	BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);

	vLower = v->getValue(i1, i2, *wLower, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All);
	vUpper = v->getValue(i1, i2, *wUpper, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All);
	vExact = v->getValue(i1, i2, (1.0f - lambda) * *wLower + lambda * *wUpper, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::U);
	BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);


	vLower = v->getValue(i1, i2, *wLower, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All);
	vUpper = v->getValue(i1, i2, *wUpper, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All);
	vExact = v->getValue(i1, i2, (1.0f - lambda) * *wLower + lambda * *wUpper, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::None);
	BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);

	vLower = v->getValue(i1, i2, *wExact, *wLower, *wExact, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All);
	vUpper = v->getValue(i1, i2, *wExact, *wUpper, *wExact, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All);
	vExact = v->getValue(i1, i2, *wExact, (1.0f - lambda) * *wLower + lambda * *wUpper, *wExact, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::S);
	BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);

	vLower = v->getValue(i1, i2, *wExact, *wExact, *wLower, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All);
	vUpper = v->getValue(i1, i2, *wExact, *wExact, *wUpper, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All);
	vExact = v->getValue(i1, i2, *wExact, *wExact, (1.0f - lambda) * *wLower + lambda * *wUpper, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::T);
	BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);

	vLower = v->getValue(i1, i2, *wLower, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All);
	vUpper = v->getValue(i1, i2, *wUpper, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All);
	vExact = v->getValue(i1, i2, (1.0f - lambda) * *wLower + lambda * *wUpper, *wExact, *wExact, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::U);
	BOOST_CHECK_CLOSE((1.0f - lambda) * vLower + lambda * vUpper, vExact, 0.0001);
}

BOOST_AUTO_TEST_CASE(getValueBuffered)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Spin) = float(i);
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Density) = float(i) + 1.0f;
	}

	SublatticeIterator i1 = FrgCommon::lattice().getBasis();
	SublatticeIterator i2 = ++FrgCommon::lattice().getRange(i1);

	//spin
	auto ab = v->generateAccessBuffer(1.1f, 2.2f, 3.3f);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Spin, ab), v->getValue(i1, i2, 1.1f, 2.2f, 3.3f, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::None), 0.0001); 
	
	ab = v->generateAccessBuffer(2.0f, 3.0f, 4.0f);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Spin, ab), v->getValue(i1, i2, 2.0f, 3.0f, 4.0f, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::All), 0.0001);

	auto ab2 = v->generateAccessBuffer(2.0f, 3.5f, 4.5f, SU2VertexTwoParticle::FrequencyChannel::S);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Spin, ab2), v->getValue(i1, i2, 2.0f, 3.5f, 4.5f, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::S), 0.0001);

	ab2 = v->generateAccessBuffer(2.5f, 3.0f, 4.5f, SU2VertexTwoParticle::FrequencyChannel::T);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Spin, ab2), v->getValue(i1, i2, 2.5f, 3.0f, 4.5f, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::T), 0.0001);

	ab2 = v->generateAccessBuffer(2.5f, 3.5f, 4.0f, SU2VertexTwoParticle::FrequencyChannel::U);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Spin, ab2), v->getValue(i1, i2, 2.5f, 3.5f, 4.0f, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::U), 0.0001);

	//density
	ab = v->generateAccessBuffer(1.1f, 2.2f, 3.3f);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Density, ab), v->getValue(i1, i2, 1.1f, 2.2f, 3.3f, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::None), 0.0001); 
	
	ab = v->generateAccessBuffer(2.0f, 3.0f, 4.0f);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Density, ab), v->getValue(i1, i2, 2.0f, 3.0f, 4.0f, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::All), 0.0001);

	ab2 = v->generateAccessBuffer(2.0f, 3.5f, 4.5f, SU2VertexTwoParticle::FrequencyChannel::S);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Density, ab2), v->getValue(i1, i2, 2.0f, 3.5f, 4.5f, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::S), 0.0001);

	ab2 = v->generateAccessBuffer(2.5f, 3.0f, 4.5f, SU2VertexTwoParticle::FrequencyChannel::T);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Density, ab2), v->getValue(i1, i2, 2.5f, 3.0f, 4.5f, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::T), 0.0001);

	ab2 = v->generateAccessBuffer(2.5f, 3.5f, 4.0f, SU2VertexTwoParticle::FrequencyChannel::U);
	BOOST_CHECK_CLOSE(v->getValue(i1, i2, SU2VertexTwoParticle::Symmetry::Density, ab2), v->getValue(i1, i2, 2.5f, 3.5f, 4.0f, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::U), 0.0001);
}

BOOST_AUTO_TEST_CASE(getValueLocal)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Spin) = float(i);
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Density) = float(i) + 1.0f;
	}

	//spin
	auto ab = v->generateAccessBuffer(1.1f, 2.2f, 3.3f);
	BOOST_CHECK_CLOSE(v->getValueLocal(SU2VertexTwoParticle::Symmetry::Spin, ab), v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), 1.1f, 2.2f, 3.3f, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::None), 0.0001);

	//density
	ab = v->generateAccessBuffer(1.1f, 2.2f, 3.3f);
	BOOST_CHECK_CLOSE(v->getValueLocal(SU2VertexTwoParticle::Symmetry::Density, ab), v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), 1.1f, 2.2f, 3.3f, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::None), 0.0001);
}

BOOST_AUTO_TEST_CASE(getValueSuperbundle)
{
	for (int i = 0; i < v->size; ++i)
	{
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Spin) = float(i);
		v->getValueRef(i, SU2VertexTwoParticle::Symmetry::Density) = float(i) + 1.0f;
	}

	auto ab = v->generateAccessBuffer(1.1f, 2.2f, 3.3f);
	ValueSuperbundle<float, 2> b(FrgCommon::lattice().size);
	v->getValueSuperbundle(ab, b);

	for (int rid = 0; rid < FrgCommon::lattice().size; ++rid)
	{
		BOOST_CHECK_CLOSE(b.bundle(0)[rid], v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().fromParametrization(rid), 1.1f, 2.2f, 3.3f, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::None), 0.0001);
		BOOST_CHECK_CLOSE(b.bundle(1)[rid], v->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().fromParametrization(rid), 1.1f, 2.2f, 3.3f, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::None), 0.0001);
	}
}

BOOST_AUTO_TEST_SUITE_END();