#define BOOST_TEST_MODULE "LoadManagerTest"
#include <chrono>
#include <thread>
#include <boost/test/included/unit_test.hpp>
#include "lib/LoadManager.hpp"

#ifndef DISABLE_MPI
#include "mpi.h"
#endif

struct MPIFixture
{
	MPIFixture()
	{
		#ifndef DISABLE_MPI
		int argc = boost::unit_test::framework::master_test_suite().argc;
		char **argv = boost::unit_test::framework::master_test_suite().argv;
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		#endif
	}

	~MPIFixture()
	{
		#ifndef DISABLE_MPI
		MPI_Finalize();
		#endif
	}

	static int rank;
};
int MPIFixture::rank = 0;

BOOST_GLOBAL_FIXTURE(MPIFixture);

struct LoadManagerFixture
{
	LoadManagerFixture()
	{
		m = HMP::newLoadManager(0);
	}

	~LoadManagerFixture()
	{
		delete m;
	}

	HMP::LoadManager *m;
};

BOOST_FIXTURE_TEST_SUITE(LoadManagerTestMPI, LoadManagerFixture);

BOOST_AUTO_TEST_CASE(MasterStackExplicit)
{
	const int dataLength = 16;
	float data1[dataLength];
	float data2[dataLength];

	auto resetdata = [&]()->void {
		for (int i = 0; i < dataLength; ++i)
		{
			data1[i] = float(i);
			data2[i] = -float(i);
		}
	};

	std::function<float(int)> calculator1 = [](int n)->float { std::this_thread::sleep_for(std::chrono::milliseconds(50)); return float(n * n); };
	std::function<float(int)> calculator2 = [](int n)->float { std::this_thread::sleep_for(std::chrono::milliseconds(50)); return float(n * n * n); };

	HMP::StackIdentifier stack1 = m->addMasterStackExplicit(&data1[0], dataLength, calculator1);
	HMP::StackIdentifier stack2 = m->addMasterStackExplicit(&data2[0], dataLength, calculator2);

	resetdata();
	m->calculate(stack1);
	m->broadcast(stack1);
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data1[i], float(i * i));
	m->calculate(stack2);
	m->broadcast(stack2);
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data2[i], float(i * i * i));

	resetdata();
	m->calculateAll();
	m->broadcastAll();
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data1[i], float(i * i));
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data2[i], float(i * i * i));
}

BOOST_AUTO_TEST_CASE(MasterStackExplicitAutobroadcast)
{
	const int dataLength = 16;
	float data1[dataLength];
	float data2[dataLength];

	auto resetdata = [&]()->void {
		for (int i = 0; i < dataLength; ++i)
		{
			data1[i] = float(i);
			data2[i] = -float(i);
		}
	};

	std::function<float(int)> calculator1 = [](int n)->float { std::this_thread::sleep_for(std::chrono::milliseconds(50)); return float(n * n); };
	std::function<float(int)> calculator2 = [](int n)->float { std::this_thread::sleep_for(std::chrono::milliseconds(50)); return float(n * n * n); };

	m->addMasterStackExplicit(&data1[0], dataLength, calculator1, 1, 8, true);
	m->addMasterStackExplicit(&data2[0], dataLength, calculator2, 1, 8, true);

	resetdata();
	m->calculateAll();
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data1[i], float(i * i));
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data2[i], float(i * i * i));
}

BOOST_AUTO_TEST_CASE(MasterStackImplicit)
{
	const int dataLength = 16;
	float data1[dataLength];

	auto resetdata = [&]()->void {
		for (int i = 0; i < dataLength; ++i)
		{
			data1[i] = float(i);
		}
	};

	std::function<void(int)> calculator1 = [&data1](int n)->void { std::this_thread::sleep_for(std::chrono::milliseconds(50)); data1[n] = float(n * n); };

	HMP::StackIdentifier stack1 = m->addMasterStackImplicit(&data1[0], dataLength, calculator1, 1, 1, 8, true);

	resetdata();
	m->calculate(stack1);
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data1[i], float(i * i));
}

BOOST_AUTO_TEST_CASE(MasterStackImplicitTyplemult)
{
	const int dataLength = 8;
	const int dataMultiplicity = 2;
	float data1[dataLength * dataMultiplicity];

	auto resetdata = [&]()->void {
		for (int i = 0; i < dataLength * dataMultiplicity; ++i)
		{
			data1[i] = float(i);
		}
	};

	std::function<void(int)> calculator1 = [&data1](int n)->void { std::this_thread::sleep_for(std::chrono::milliseconds(50)); data1[2*n] = float(4 * n * n); data1[2*n+1] = float((2 * n + 1) * (2 * n + 1)); };

	HMP::StackIdentifier stack1 = m->addMasterStackImplicit(&data1[0], dataLength, calculator1, dataMultiplicity, 1, 4, true);

	resetdata();
	m->calculate(stack1);
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data1[i], float(i * i));
}

BOOST_AUTO_TEST_CASE(SlaveStack)
{
	const int dataLength = 8;
	const int dataMultiplicity = 2;
	float data1[dataLength * dataMultiplicity];
	float data2[dataLength * dataMultiplicity];

	auto resetdata = [&]()->void {
		for (int i = 0; i < dataLength * dataMultiplicity; ++i)
		{
			data1[i] = float(i);
			data2[i] = float(i);
		}
	};

	std::function<void(int)> calculator1 = [&data1,&data2](int n)->void { 
		std::this_thread::sleep_for(std::chrono::milliseconds(50)); 
		data1[2 * n] = float((2 * n) * (2 * n)); 
		data1[2 * n + 1] = float((2 * n + 1) * (2 * n + 1)); 
		data2[2 * n] = float((2 * n) * (2 * n) * (2 * n));
		data2[2 * n + 1] = float((2 * n + 1) * (2 * n + 1) * (2 * n + 1));
	};

	HMP::StackIdentifier stack1 = m->addMasterStackImplicit(&data1[0], dataLength, calculator1, dataMultiplicity, 1, 4, true);
	m->addSlaveStack(&data2[0], dataLength, stack1, dataMultiplicity);

	BOOST_CHECK_THROW(m->addSlaveStack(&data2[0], dataLength, stack1, 2 * dataMultiplicity), Exception);

	resetdata();
	m->calculate(stack1);
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data1[i], float(i * i));
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data2[i], float(i * i * i));
}

BOOST_AUTO_TEST_CASE(PassiveStack)
{
	const int dataLength = 16;
	float data1[dataLength];

	auto resetdata = [&]()->void {
		for (int i = 0; i < dataLength; ++i)
		{
			data1[i] = float(i);
		}
	};

	m->addPassiveStack(&data1[0], dataLength);

	resetdata();
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data1[i], float(i));

	if (MPIFixture::rank == 0)
	{
		for (int i = 0; i < dataLength; ++i) data1[i] = float(i * i);
	}
	m->broadcastAll();
	for (int i = 0; i < dataLength; ++i) BOOST_CHECK_EQUAL(data1[i], float(i * i));
}

BOOST_AUTO_TEST_SUITE_END();