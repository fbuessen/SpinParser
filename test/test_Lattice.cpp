#define BOOST_TEST_MODULE "LatticeTest"
#include <boost/test/included/unit_test.hpp>
#include "lib/Log.hpp"
#include "LatticeModelFactory.hpp"

typedef std::tuple<int, int, int, int> SiteParameters;
BOOST_TEST_DONT_PRINT_LOG_VALUE(SiteParameters);
BOOST_TEST_DONT_PRINT_LOG_VALUE(geometry::Vec3<double>);
BOOST_TEST_DONT_PRINT_LOG_VALUE(SpinComponent);

struct SquareLatticeFixture
{
	SquareLatticeFixture()
	{
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
		l = product.first;
		delete product.second;
	}

	~SquareLatticeFixture()
	{
		delete l;
	}

	Lattice *l;
};

struct HoneycombLatticeFixture
{
	HoneycombLatticeFixture()
	{
		LatticeModelFactory::LatticeUnitCell uc;
		uc.basisSites.push_back(geometry::Vec3<double>(0.0, 0.0, 0.0));
		uc.basisSites.push_back(geometry::Vec3<double>(1.0, 0.0, 0.0));
		uc.latticeVectors.push_back(geometry::Vec3<double>(1.5, 0.5 * sqrt(3), 0.0));
		uc.latticeVectors.push_back(geometry::Vec3<double>(1.5, -0.5 * sqrt(3), 0.0));
		uc.latticeVectors.push_back(geometry::Vec3<double>(0.0, 0.0, 1.0));
		uc.latticeBonds.push_back(LatticeModelFactory::LatticeBond(0, 1, 0, 0, 0));
		uc.latticeBonds.push_back(LatticeModelFactory::LatticeBond(1, 0, 1, 0, 0));
		uc.latticeBonds.push_back(LatticeModelFactory::LatticeBond(1, 0, 0, 1, 0));

		LatticeModelFactory::SpinModelUnitCell model;
		LatticeModelFactory::SpinInteraction i1(LatticeModelFactory::LatticeSite(0, 0, 0, 0), LatticeModelFactory::LatticeSite(0, 0, 0, 1));
		i1.interactionStrength[0][0] = 0.1f;
		model.interactions.push_back(i1);
		LatticeModelFactory::SpinInteraction i2(LatticeModelFactory::LatticeSite(0, 0, 0, 0), LatticeModelFactory::LatticeSite(0, -1, 0, 1));
		i2.interactionStrength[1][1] = 0.2f;
		model.interactions.push_back(i2);
		LatticeModelFactory::SpinInteraction i3(LatticeModelFactory::LatticeSite(0, 0, 0, 0), LatticeModelFactory::LatticeSite(-1, 0, 0, 1));
		i3.interactionStrength[2][2] = 0.2f;
		model.interactions.push_back(i3);

		Log::log << Log::setDisplayLogLevel(Log::LogLevel::None);
		std::pair<Lattice *, SpinModel *> product = LatticeModelFactory::newLatticeModel(uc, model, 3);
		l = product.first;
		delete product.second;
	}

	~HoneycombLatticeFixture()
	{
		delete l;
	}

	Lattice *l;
};

BOOST_AUTO_TEST_SUITE(LatticeTest);

BOOST_FIXTURE_TEST_CASE(SquareLatticeIterate, SquareLatticeFixture)
{
	BOOST_CHECK_EQUAL((l->getSiteParameters(l->getBasis())), (SiteParameters(0, 0, 0, 0)));
	BOOST_CHECK_EQUAL((l->getSiteParameters(l->zero())), (SiteParameters(0, 0, 0, 0)));

	int n = 0;
	for (auto i = l->getRange(l->getBasis()); i != l->end(); ++i) ++n;
	BOOST_CHECK_EQUAL(n, 25);

	BOOST_CHECK_EQUAL(l->size, 6);
};

BOOST_FIXTURE_TEST_CASE(SquareLatticeGeometry, SquareLatticeFixture)
{
	for (auto i = l->getRange(l->getBasis()); i != l->end(); ++i)
	{
		if (l->getSiteParameters(i) == SiteParameters(2, 1, 0, 0)) BOOST_CHECK_EQUAL(l->getSitePosition(i), geometry::Vec3<double>(2.0, 1.0, 0.0));
		if (l->getSiteParameters(i) == SiteParameters(-2, 1, 0, 0)) BOOST_CHECK_EQUAL(l->getSitePosition(i), geometry::Vec3<double>(-2.0, 1.0, 0.0));
	}
};

BOOST_FIXTURE_TEST_CASE(SquareLatticeSymmetry, SquareLatticeFixture)
{
	std::vector<SiteParameters> i1List = { SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0) };
	std::vector<SiteParameters> i2List = { SiteParameters(0,0,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-2,0,0),SiteParameters(-1,-1,0,0),SiteParameters(0,-3,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,0,0,0),SiteParameters(1,1,0,0),SiteParameters(2,0,0,0),SiteParameters(1,-1,0,0),SiteParameters(-1,1,0,0),SiteParameters(-2,0,0,0),SiteParameters(0,2,0,0),SiteParameters(0,1,0,0),SiteParameters(1,0,0,0),SiteParameters(2,1,0,0),SiteParameters(2,-1,0,0),SiteParameters(1,2,0,0),SiteParameters(3,0,0,0),SiteParameters(1,-2,0,0),SiteParameters(-2,1,0,0),SiteParameters(-2,-1,0,0),SiteParameters(-1,2,0,0),SiteParameters(-3,0,0,0),SiteParameters(0,3,0,0),SiteParameters(0,0,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-2,0,0),SiteParameters(-1,-1,0,0),SiteParameters(-1,0,0,0),SiteParameters(1,1,0,0),SiteParameters(2,0,0,0),SiteParameters(1,-1,0,0),SiteParameters(-1,1,0,0),SiteParameters(-2,0,0,0),SiteParameters(0,2,0,0),SiteParameters(0,1,0,0),SiteParameters(1,0,0,0),SiteParameters(2,1,0,0),SiteParameters(2,-1,0,0),SiteParameters(1,2,0,0),SiteParameters(3,0,0,0),SiteParameters(1,-2,0,0) };
	std::vector<SiteParameters> targetList = { SiteParameters(0,0,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-2,0,0),SiteParameters(-1,-1,0,0),SiteParameters(0,-3,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-1,0,0),SiteParameters(-1,-1,0,0),SiteParameters(0,-2,0,0),SiteParameters(-1,-1,0,0),SiteParameters(-1,-1,0,0),SiteParameters(0,-2,0,0),SiteParameters(0,-2,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-1,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-3,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-3,0,0),SiteParameters(0,-3,0,0),SiteParameters(0,-1,0,0),SiteParameters(-1,-1,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-2,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-1,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-3,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-1,0,0),SiteParameters(0,0,0,0),SiteParameters(-1,-1,0,0),SiteParameters(-1,-1,0,0),SiteParameters(0,-2,0,0),SiteParameters(0,-2,0,0),SiteParameters(0,-2,0,0) };

	auto fromParameters = [&](const SiteParameters &p)->LatticeIterator
	{
		for (auto i = l->begin(); i != l->end(); ++i)
		{
			if (l->getSiteParameters(i) == p) return i;
		}
		return l->end();
	};

	for (int n = 0; n < int(i1List.size()); ++n)
	{
		LatticeIterator i1 = fromParameters(i1List[n]);
		LatticeIterator i2 = fromParameters(i2List[n]);
		LatticeIterator target = fromParameters(targetList[n]);

		BOOST_CHECK_EQUAL(l->symmetryTransform(i1, i2), l->symmetryTransform(l->zero(), target));
	}
};

BOOST_FIXTURE_TEST_CASE(SquareLatticeOverlap, SquareLatticeFixture)
{
	std::vector<SiteParameters> i1List = { SiteParameters(0,0,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-2,0,0),SiteParameters(-1,-1,0,0),SiteParameters(0,-3,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-1,0,0),SiteParameters(-1,-1,0,0),SiteParameters(-1,-1,0,0),SiteParameters(0,-2,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-1,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-3,0,0) };
	std::vector<SiteParameters> i2List = { SiteParameters(-1,-1,0,0),SiteParameters(0,-1,0,0),SiteParameters(-1,-1,0,0),SiteParameters(0,0,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-2,0,0),SiteParameters(0,-2,0,0),SiteParameters(-1,-1,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-3,0,0),SiteParameters(-1,-2,0,0),SiteParameters(-1,-2,0,0),SiteParameters(0,-1,0,0),SiteParameters(0,-3,0,0),SiteParameters(-1,-2,0,0) };

	auto fromParameters = [&](const SiteParameters &p)->LatticeIterator
	{
		for (auto i = l->begin(); i != l->end(); ++i)
		{
			if (l->getSiteParameters(i) == p) return i;
		}
		return l->end();
	};
	auto findPairOnce = [&](int rid1, int rid2)->bool
	{
		for (int n = 0; n < int(i1List.size()); ++n)
		{
			LatticeIterator i1 = fromParameters(i1List[n]);
			LatticeIterator i2 = fromParameters(i2List[n]);
			if (rid1 == l->symmetryTransform(l->zero(), i1) && rid2 == l->symmetryTransform(l->zero(), i2))
			{
				i1List.erase(i1List.begin() + n);
				i2List.erase(i2List.begin() + n);
				return true;
			}
		}
		return false;
	};

	const LatticeOverlap &o = l->getOverlap(l->symmetryTransform(l->zero(), fromParameters(SiteParameters(1, 1, 0, 0))));
	BOOST_CHECK_EQUAL(o.size, 18);

	for (int n = 0; n < o.size; ++n)
	{
		BOOST_CHECK(findPairOnce(o.rid1[n], o.rid2[n]));
	}
};

BOOST_FIXTURE_TEST_CASE(HoneycombLatticeIterate, HoneycombLatticeFixture)
{
	BOOST_CHECK_EQUAL((l->getSiteParameters(l->getBasis())), (SiteParameters(0, 0, 0, 0)));
	BOOST_CHECK_EQUAL((l->getSiteParameters(++l->getBasis())), (SiteParameters(0, 0, 0, 1)));
	BOOST_CHECK_EQUAL((l->getSiteParameters(l->zero())), (SiteParameters(0, 0, 0, 0)));

	int n = 0;
	for (auto i = l->getRange(l->getBasis()); i != l->end(); ++i) ++n;
	BOOST_CHECK_EQUAL(n, 19);

	BOOST_CHECK_EQUAL(l->size, 11);
};

BOOST_FIXTURE_TEST_CASE(HoneycombLatticeGeometry, HoneycombLatticeFixture)
{
	for (auto i = l->getRange(l->getBasis()); i != l->end(); ++i)
	{
		if (l->getSiteParameters(i) == SiteParameters(1, 0, 0, 0)) BOOST_CHECK_EQUAL(l->getSitePosition(i), geometry::Vec3<double>(1.5, 0.5 * sqrt(3), 0.0));
		else if (l->getSiteParameters(i) == SiteParameters(-1, 0, 0, 0)) BOOST_CHECK_EQUAL(l->getSitePosition(i), geometry::Vec3<double>(-1.5, -0.5 * sqrt(3), 0.0));
		else if (l->getSiteParameters(i) == SiteParameters(1, 0, 0, 1)) BOOST_CHECK_EQUAL(l->getSitePosition(i), geometry::Vec3<double>(2.5, 0.5 * sqrt(3), 0.0));
		else if (l->getSiteParameters(i) == SiteParameters(-1, 0, 0, 1)) BOOST_CHECK_EQUAL(l->getSitePosition(i), geometry::Vec3<double>(-0.5, -0.5 * sqrt(3), 0.0));
	}
};

BOOST_FIXTURE_TEST_CASE(HoneycombLatticeSymmetry, HoneycombLatticeFixture)
{
	std::vector<SiteParameters> i1List = { SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(0,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0),SiteParameters(1,0,0,0) };
	std::vector<SiteParameters> i2List = { SiteParameters(0,0,0,0),SiteParameters(0,0,0,1),SiteParameters(-1,0,0,1),SiteParameters(1,0,0,0),SiteParameters(-1,0,0,0),SiteParameters(-1,1,0,0),SiteParameters(1,0,0,1),SiteParameters(1,-1,0,1),SiteParameters(-2,0,0,1),SiteParameters(-1,-1,0,1),SiteParameters(-2,1,0,1),SiteParameters(0,1,0,0),SiteParameters(0,1,0,1),SiteParameters(-1,1,0,1),SiteParameters(0,-1,0,0),SiteParameters(1,-1,0,0),SiteParameters(0,-1,0,1),SiteParameters(0,-2,0,1),SiteParameters(1,-2,0,1),SiteParameters(0,0,0,0),SiteParameters(0,0,0,1),SiteParameters(-1,0,0,1),SiteParameters(1,0,0,0),SiteParameters(1,0,0,1),SiteParameters(1,-1,0,1),SiteParameters(0,1,0,0),SiteParameters(0,1,0,1),SiteParameters(-1,1,0,1),SiteParameters(1,-1,0,0),SiteParameters(0,-1,0,1),SiteParameters(1,-2,0,1) };
	std::vector<SiteParameters> targetList = { SiteParameters(0,0,0,0),SiteParameters(0,0,0,1),SiteParameters(-1,0,0,1),SiteParameters(1,0,0,0),SiteParameters(-1,0,0,0),SiteParameters(-1,1,0,0),SiteParameters(1,0,0,1),SiteParameters(1,-1,0,1),SiteParameters(-2,0,0,1),SiteParameters(-1,-1,0,1),SiteParameters(-2,1,0,1),SiteParameters(1,0,0,0),SiteParameters(1,0,0,1),SiteParameters(1,-1,0,1),SiteParameters(-1,0,0,0),SiteParameters(-1,1,0,0),SiteParameters(-1,0,0,1),SiteParameters(-2,0,0,1),SiteParameters(-2,1,0,1),SiteParameters(-1,0,0,0),SiteParameters(-1,0,0,1),SiteParameters(-2,0,0,1),SiteParameters(0,0,0,0),SiteParameters(0,0,0,1),SiteParameters(-1,0,0,1),SiteParameters(-1,1,0,0),SiteParameters(1,-1,0,1),SiteParameters(-2,1,0,1),SiteParameters(-1,0,0,0),SiteParameters(-1,-1,0,1),SiteParameters(-2,0,0,1) };
	std::vector<SpinComponent> sxTarget = { SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X };
	std::vector<SpinComponent> syTarget = { SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Z,SpinComponent::Y,SpinComponent::Z,SpinComponent::Y,SpinComponent::Z,SpinComponent::Y,SpinComponent::Z };
	std::vector<SpinComponent> szTarget = { SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Y,SpinComponent::Z,SpinComponent::Y,SpinComponent::Z,SpinComponent::Y,SpinComponent::Z,SpinComponent::Y };

	auto fromParameters = [&](const SiteParameters &p)->LatticeIterator
	{
		for (auto i = l->begin(); i != l->end(); ++i)
		{
			if (l->getSiteParameters(i) == p) return i;
		}
		return l->end();
	};

	for (int n = 0; n < int(i1List.size()); ++n)
	{
		LatticeIterator i1 = fromParameters(i1List[n]);
		LatticeIterator i2 = fromParameters(i2List[n]);
		LatticeIterator target = fromParameters(targetList[n]);

		SpinComponent s1 = SpinComponent::X;
		SpinComponent s2 = SpinComponent::Y;
		SpinComponent s3 = SpinComponent::Z;
		int t = l->symmetryTransform(i1, i2, s1, s2, s3);

		SpinComponent s1p = sxTarget[n];
		SpinComponent s2p = syTarget[n];
		SpinComponent s3p = szTarget[n];
		int tp = l->symmetryTransform(l->zero(), target, s1p, s2p, s3p);

		BOOST_CHECK_EQUAL(t, tp);
		BOOST_CHECK_EQUAL(s1, s1p);
		BOOST_CHECK_EQUAL(s2, s2p);
		BOOST_CHECK_EQUAL(s3, s3p);
	}
};

BOOST_FIXTURE_TEST_CASE(HoneycombLatticeOverlap, HoneycombLatticeFixture)
{
	std::vector<SiteParameters> i1List = { SiteParameters(0,0,0,0),SiteParameters(0,0,0,1),SiteParameters(1,0,0,0),SiteParameters(1,0,0,1),SiteParameters(1,-1,0,1),SiteParameters(1,0,0,0),SiteParameters(1,0,0,1),SiteParameters(-1,1,0,0) };
	std::vector<SiteParameters> i2List = { SiteParameters(1,0,0,1),SiteParameters(-1,0,0,0),SiteParameters(0,0,0,1),SiteParameters(0,0,0,0),SiteParameters(-1,0,0,0),SiteParameters(1,-1,0,1),SiteParameters(-1,1,0,0),SiteParameters(1,0,0,1) };
	std::vector<SpinComponent> s1xList = { SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X };
	std::vector<SpinComponent> s1yList = { SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z };
	std::vector<SpinComponent> s1zList = { SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y };
	std::vector<SpinComponent> s2xList = { SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X,SpinComponent::X };
	std::vector<SpinComponent> s2yList = { SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Y,SpinComponent::Z,SpinComponent::Y,SpinComponent::Y,SpinComponent::Z };
	std::vector<SpinComponent> s2zList = { SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Z,SpinComponent::Y,SpinComponent::Z,SpinComponent::Z,SpinComponent::Y };

	auto fromParameters = [&](const SiteParameters &p)->LatticeIterator
	{
		for (auto i = l->begin(); i != l->end(); ++i)
		{
			if (l->getSiteParameters(i) == p) return i;
		}
		return l->end();
	};
	auto findPairOnce = [&](int rid1, int rid2, SpinComponent tx1, SpinComponent ty1, SpinComponent tz1, SpinComponent tx2, SpinComponent ty2, SpinComponent tz2)->bool
	{
		for (int n = 0; n < int(i1List.size()); ++n)
		{
			LatticeIterator i1 = fromParameters(i1List[n]);
			LatticeIterator i2 = fromParameters(i2List[n]);

			SpinComponent sx1 = s1xList[n];
			SpinComponent sy1 = s1yList[n];
			SpinComponent sz1 = s1zList[n];
			SpinComponent sx2 = s2xList[n];
			SpinComponent sy2 = s2yList[n];
			SpinComponent sz2 = s2zList[n];
			int i1p = l->symmetryTransform(l->zero(), i1, sx1, sy1, sz1);
			int i2p = l->symmetryTransform(l->zero(), i2, sx2, sy2, sz2);

			if (rid1 == i1p && rid2 == i2p && tx1 == sx1 && ty1 == sy1 && tz1 == sz1 && tx2 == sx2 && ty2 == sy2 && tz2 == sz2)
			{
				i1List.erase(i1List.begin() + n);
				i2List.erase(i2List.begin() + n);
				s1xList.erase(s1xList.begin() + n);
				s1yList.erase(s1yList.begin() + n);
				s1zList.erase(s1zList.begin() + n);
				s2xList.erase(s2xList.begin() + n);
				s2yList.erase(s2yList.begin() + n);
				s2zList.erase(s2zList.begin() + n);
				return true;
			}
		}
		return false;
	};

	const LatticeOverlap &o = l->getOverlap(l->symmetryTransform(l->zero(), fromParameters(SiteParameters(1, 0, 0, 1))));
	BOOST_CHECK_EQUAL(o.size, 8);

	for (int n = 0; n < o.size; ++n)
	{
		BOOST_CHECK(findPairOnce(o.rid1[n], o.rid2[n], o.transformedX1[n], o.transformedY1[n], o.transformedZ1[n], o.transformedX2[n], o.transformedY2[n], o.transformedZ2[n]));
	}
};

BOOST_AUTO_TEST_SUITE_END();