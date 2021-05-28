/**
 * @file SU2FrgCore.cpp
 * @author Finn Lasse Buessen
 * @brief FrgCore implementation for SU(2) models.
 * 
 * @copyright Copyright (c) 2020
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include "lib/InputParser.hpp"
#include "lib/Integrator.hpp"
#include "SpinParser.hpp"
#include "SU2FrgCore.hpp"
#include "SU2EffectiveAction.hpp"

SU2FrgCore::SU2FrgCore(const SpinModel &spinModel, const std::vector<Measurement *> &measurements, const std::map<std::string, std::string> &options) : FrgCore(measurements)
{
	//init options
	spinLength = 0.5;
	normalization = NAN;

	for (auto option : options)
	{
		if (option.first == "spin") spinLength = InputParser::stringToFloat(option.second);
		else if (option.first == "normalization") normalization = InputParser::stringToFloat(option.second);
		else throw Exception(Exception::Type::InitializationError, "Unknown spin model option '" + option.first + "'.");
	}
	if (std::isnan(normalization)) normalization = 2.0f * spinLength;

	Log::log << Log::LogLevel::Info << "FRG core spin length S is set to " << spinLength << "." << Log::endl;
	Log::log << Log::LogLevel::Info << "FRG core energy normalization is set to " << normalization << "." << Log::endl;

	//init data
	_flowingFunctional = new SU2EffectiveAction(*FrgCommon::cutoff().begin(), spinModel, this);
	_flow = new SU2EffectiveAction();

	//init loadManager
	//stack0
	dataStacks[0] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		&_flowingFunctional->cutoff,
		1);
	//stack1
	dataStacks[1] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexSingleParticle->_data,
		static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexSingleParticle->size);
	//stack2
	dataStacks[2] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataDD,
		static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size);
	//stack3
	dataStacks[3] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataSS,
		static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size);
	//stack4
	dataStacks[4] = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		&_flow->cutoff,
		1,
		[&](int) { _flow->cutoff = static_cast<SU2EffectiveAction *>(_flowingFunctional)->cutoff; });
	//stack5
	dataStacks[5] = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		static_cast<SU2EffectiveAction *>(_flow)->vertexSingleParticle->_data,
		static_cast<SU2EffectiveAction *>(_flow)->vertexSingleParticle->size,
		[&](int x) { _calculateVertexSingleParticle(x); },
		1,
		1,
		1);
	//stack6
	dataStacks[6] = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		static_cast<SU2EffectiveAction *>(_flow)->vertexTwoParticle->_dataDD,
		static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->sizeFrequency,
		[&](int x) { _calculateVertexTwoParticle(x); },
		FrgCommon::lattice().size,
		FrgCommon::frequency().size);
	//stack7
	dataStacks[7] = SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		static_cast<SU2EffectiveAction *>(_flow)->vertexTwoParticle->_dataSS,
		static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->sizeFrequency,
		dataStacks[6],
		FrgCommon::lattice().size);
}

SU2FrgCore::~SU2FrgCore()
{
	delete _flowingFunctional;
	delete _flow;
}

void SU2FrgCore::computeStep()
{
	//update cutoff and broadcast
	SpinParser::spinParser()->getLoadManager()->calculate(dataStacks[4]);
	SpinParser::spinParser()->getLoadManager()->broadcast(dataStacks[4]);
	//calculate 1-particle vertices and broadcast (required for Katanin calculation)
	SpinParser::spinParser()->getLoadManager()->calculate(dataStacks[5]);
	SpinParser::spinParser()->getLoadManager()->broadcast(dataStacks[5]);
	//calculate 2-particle vertices and managed measurements
	std::vector<int> managedMeasurementStacks;
	for (auto m = _measurements.begin(); m != _measurements.end(); ++m)
	{
		if ((*m)->isLoadManaged())
		{
			auto s = (*m)->getLoadManagedStacks();
			managedMeasurementStacks.insert(managedMeasurementStacks.end(), s.begin(), s.end());
		}
	}
	managedMeasurementStacks.push_back(dataStacks[6]);
	SpinParser::spinParser()->getLoadManager()->calculate(managedMeasurementStacks.data(), int(managedMeasurementStacks.size()));
}

void SU2FrgCore::finalizeStep(float newCutoff)
{
	//determine cutoff set
	float cutoffStep = newCutoff - _flowingFunctional->cutoff;

	//set new cutoff value
	_flowingFunctional->cutoff = newCutoff;

	//add _flow to single particle vertex
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexSingleParticle->size; ++i) static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexSingleParticle->_data[i] += cutoffStep * static_cast<SU2EffectiveAction *>(_flow)->vertexSingleParticle->_data[i];

	//add _flow to two particle vertex
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size; ++i) static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataDD[i] += cutoffStep * static_cast<SU2EffectiveAction *>(_flow)->vertexTwoParticle->_dataDD[i];
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size; ++i) static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataSS[i] += cutoffStep * static_cast<SU2EffectiveAction *>(_flow)->vertexTwoParticle->_dataSS[i];

	//broadcast updated effective action
	SpinParser::spinParser()->getLoadManager()->broadcast({ dataStacks[0], dataStacks[1], dataStacks[2], dataStacks[3] });
}

void SU2FrgCore::_calculateVertexSingleParticle(const int iterator)
{
	float cutoff = _flowingFunctional->cutoff;
	SU2VertexSingleParticle *v2 = static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexSingleParticle;
	SU2VertexTwoParticle *v4 = static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle;

	float w;
	v2->expandIterator(iterator, w);
	float v2CurrentValue = 0;

	//term1
	float sum = 0;
	for (auto j = FrgCommon::lattice().getRange(0); j != FrgCommon::lattice().end(); ++j)
	{
		sum += v4->getValue(FrgCommon::lattice().zero(), j, w + cutoff, 0.0f, w - cutoff, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::None);
		sum -= v4->getValue(FrgCommon::lattice().zero(), j, w - cutoff, 0.0f, w + cutoff, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::None);
	}
	//additional prefactor 2.0 * spinLength from generalization to arbitrary spin
	v2CurrentValue -= 4.0f * spinLength * sum;

	//term2
	v2CurrentValue += 0.75f * (v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SU2VertexTwoParticle::Symmetry::Spin, SU2VertexTwoParticle::FrequencyChannel::None));

	//term3
	v2CurrentValue += (v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SU2VertexTwoParticle::Symmetry::Density, SU2VertexTwoParticle::FrequencyChannel::None));

	//prefactor
	v2CurrentValue /= (2.0f * float(M_PI) * (cutoff + v2->getValue(cutoff)));

	static_cast<SU2EffectiveAction *>(_flow)->vertexSingleParticle->_data[iterator] = v2CurrentValue;
}

void SU2FrgCore::_calculateVertexTwoParticle(const int iterator)
{
	float cutoff = _flowingFunctional->cutoff;
	SU2VertexSingleParticle *v2 = static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexSingleParticle;
	SU2VertexTwoParticle *v4 = static_cast<SU2EffectiveAction *>(_flowingFunctional)->vertexTwoParticle;

	float s, t, u;
	v4->expandIterator(iterator, s, t, u);

	//vertex buffers
	ValueSuperbundle<float, 2> buffer1(FrgCommon::lattice().size);
	ValueSuperbundle<float, 2> buffer2(FrgCommon::lattice().size);
	ValueSuperbundle<float, 2> bufferRPA(FrgCommon::lattice().size);
	ValueSuperbundle<float, 2> v4CurrentValue(FrgCommon::lattice().size);
	ValueSuperbundle<float, 2> stackBuffers[4] = {
		ValueSuperbundle<float, 2>(FrgCommon::lattice().size),
		ValueSuperbundle<float, 2>(FrgCommon::lattice().size),
		ValueSuperbundle<float, 2>(FrgCommon::lattice().size),
		ValueSuperbundle<float, 2>(FrgCommon::lattice().size)
	};

	//transfer frequencies
	float w1p = 0.5f * (s + t + u);
	float w1 = 0.5f * (s - t + u);
	float w2p = 0.5f * (s - t - u);
	float w2 = 0.5f * (s + t - u);

	//integrand of the ferquency integral
	auto integralKernelS = [&](const float wp, ValueSuperbundle<float, 2> &returnBuffer) -> void
	{
		//pp-ladder A and B (positive sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab0 = v4->generateAccessBuffer(s, -w1 - wp, -w2 - wp, SU2VertexTwoParticle::FrequencyChannel::S);
		const SU2VertexTwoParticleAccessBuffer<4> ab1 = v4->generateAccessBuffer(s, w1p + wp, -w2p - wp, SU2VertexTwoParticle::FrequencyChannel::S);
		//pp-ladder A and B (positive sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab2 = v4->generateAccessBuffer(s, w2 + wp, w1 + wp, SU2VertexTwoParticle::FrequencyChannel::S);
		const SU2VertexTwoParticleAccessBuffer<4> ab3 = v4->generateAccessBuffer(s, -w2p - wp, w1p + wp, SU2VertexTwoParticle::FrequencyChannel::S);

		v4->getValueSuperbundle(ab0, stackBuffers[0]);
		v4->getValueSuperbundle(ab1, stackBuffers[1]);
		v4->getValueSuperbundle(ab2, stackBuffers[2]);
		v4->getValueSuperbundle(ab3, stackBuffers[3]);

		//calculate _flow
		returnBuffer.reset();

		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multSub(0.5f, stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multSub(0.5f, stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));

		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multAdd(3.0f / 16.0f, stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multAdd(3.0f / 16.0f, stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multAdd(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multAdd(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
	};

	auto integralKernelT = [&](const float wp, ValueSuperbundle<float, 2> &returnBuffer) -> void
	{
		//RPA diagram A and B equal chalice diagram A and inverse chalice diagram B, respectively (negative sign)
		//chalice diagram A (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab0 = v4->generateAccessBuffer(w1 - wp, t, w1p + wp, SU2VertexTwoParticle::FrequencyChannel::T);
		//inverse chalice diagram B (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab1 = v4->generateAccessBuffer(w2p - wp, t, -w2 - wp, SU2VertexTwoParticle::FrequencyChannel::T);
		//chalice diagram A (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab2 = v4->generateAccessBuffer(w1p + wp, t, w1 - wp, SU2VertexTwoParticle::FrequencyChannel::T);
		//inverse chalice diagram B (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab3 = v4->generateAccessBuffer(w2 + wp, t, wp - w2p, SU2VertexTwoParticle::FrequencyChannel::T);

		v4->getValueSuperbundle(ab0, stackBuffers[0]);
		v4->getValueSuperbundle(ab1, stackBuffers[1]);
		v4->getValueSuperbundle(ab2, stackBuffers[2]);
		v4->getValueSuperbundle(ab3, stackBuffers[3]);

		//calculate _flow
		returnBuffer.reset();

		bufferRPA.reset();
		for (int rid = 0; rid < FrgCommon::lattice().size; ++rid)
		{
			//lattice bubble
			const LatticeOverlap &overlap = FrgCommon::lattice().getOverlap(rid);

			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin))[rid] += stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin))[overlap.rid1[i]] * stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin))[overlap.rid2[i]];
			}

			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density))[rid] += stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density))[overlap.rid1[i]] * stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density))[overlap.rid2[i]];
			}
		}
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(2.0f * spinLength, bufferRPA.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multAdd(8.0f * spinLength, bufferRPA.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));

		//chalice diagram B (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab4 = v4->generateAccessBuffer(w2p - wp, -w2 - wp, t, SU2VertexTwoParticle::FrequencyChannel::U);
		//inverse chalice diagram A (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab5 = v4->generateAccessBuffer(w1 - wp, -w1p - wp, -t, SU2VertexTwoParticle::FrequencyChannel::U);
		//chalice diagram B (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab6 = v4->generateAccessBuffer(w2 + wp, wp - w2p, t, SU2VertexTwoParticle::FrequencyChannel::U);
		//inverse chalice diagram A (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab7 = v4->generateAccessBuffer(w1p + wp, wp - w1, -t, SU2VertexTwoParticle::FrequencyChannel::U);

		const float valCbs = v4->getValueLocal(SU2VertexTwoParticle::Symmetry::Spin, ab4);
		const float valCbd = v4->getValueLocal(SU2VertexTwoParticle::Symmetry::Density, ab4);
		const float valICas = v4->getValueLocal(SU2VertexTwoParticle::Symmetry::Spin, ab5);
		const float valICad = v4->getValueLocal(SU2VertexTwoParticle::Symmetry::Density, ab5);
		const float valCbs2 = v4->getValueLocal(SU2VertexTwoParticle::Symmetry::Spin, ab6);
		const float valCbd2 = v4->getValueLocal(SU2VertexTwoParticle::Symmetry::Density, ab6);
		const float valICas2 = v4->getValueLocal(SU2VertexTwoParticle::Symmetry::Spin, ab7);
		const float valICad2 = v4->getValueLocal(SU2VertexTwoParticle::Symmetry::Density, ab7);

		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multSub(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), valCbd);
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), 0.25f * valCbs);
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multSub(valICad, stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(0.25f * valICas, stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multSub(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), valCbd2);
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), 0.25f * valCbs2);
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multSub(valICad2, stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(0.25f * valICas2, stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));

		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multSub(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), valCbd);
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multSub(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), 0.75f * valCbs);
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multSub(valICad, stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multSub(0.75f * valICas, stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multSub(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), valCbd2);
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multSub(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), 0.75f * valCbs2);
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multSub(valICad2, stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multSub(0.75f * valICas2, stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
	};

	auto integralKernelU = [&](const float wp, ValueSuperbundle<float, 2> &returnBuffer) -> void
	{
		//u-Channel, to be combined with P(wp, u + wp) + P(u + wp, wp)
		//ph-ladder A and B, respectively (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab0 = v4->generateAccessBuffer(w1 + wp, wp - w2p, u, SU2VertexTwoParticle::FrequencyChannel::U);
		const SU2VertexTwoParticleAccessBuffer<4> ab1 = v4->generateAccessBuffer(w1p + wp, w2 - wp, u, SU2VertexTwoParticle::FrequencyChannel::U);
		//ph-ladder A and B, respectively (negative sign)
		const SU2VertexTwoParticleAccessBuffer<4> ab2 = v4->generateAccessBuffer(w2p - wp, -w1 - wp, u, SU2VertexTwoParticle::FrequencyChannel::U);
		const SU2VertexTwoParticleAccessBuffer<4> ab3 = v4->generateAccessBuffer(w2 - wp, w1p + wp, u, SU2VertexTwoParticle::FrequencyChannel::U);

		v4->getValueSuperbundle(ab0, stackBuffers[0]);
		v4->getValueSuperbundle(ab1, stackBuffers[1]);
		v4->getValueSuperbundle(ab2, stackBuffers[2]);
		v4->getValueSuperbundle(ab3, stackBuffers[3]);

		//calculate _flow
		returnBuffer.reset();

		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(0.5f, stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(0.5f, stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)).multAdd(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));

		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multAdd(3.0f / 16.0f, stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multAdd(3.0f / 16.0f, stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multAdd(stackBuffers[0].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), stackBuffers[1].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
		returnBuffer.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)).multAdd(stackBuffers[2].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)), stackBuffers[3].bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density)));
	};

	//propagator bubble
	auto p = [&](const float w1, const float w2) -> float
	{
		return 1.0f / ((w1 + v2->getValue(w1)) *(w2 + v2->getValue(w2)));
	};

	//Katanin contribution, call only once the _flow has been fully calculated and broadcasted
	auto pKataninContribution = [&](const float w1, const float w2) ->float
	{
		float denomW1 = w1 + v2->getValue(w1);
		return static_cast<SU2EffectiveAction *>(_flow)->vertexSingleParticle->getValue(w1) / (denomW1 * denomW1 * (w2 + v2->getValue(w2)));
	};

	//begin calculation of vertices here
	//conventional contribution
	integralKernelS(cutoff, buffer1);
	v4CurrentValue.multAdd(p(cutoff, cutoff + s), buffer1);
	if (s > 2.0f * cutoff)
	{
		integralKernelS(-cutoff, buffer1);
		v4CurrentValue.multAdd(p(cutoff, cutoff - s), buffer1);
	}
	integralKernelT(cutoff, buffer1);
	v4CurrentValue.multAdd(p(cutoff, cutoff + t), buffer1);
	if (t > 2.0f * cutoff)
	{
		integralKernelT(-cutoff, buffer1);
		v4CurrentValue.multAdd(p(cutoff, cutoff - t), buffer1);
	}
	integralKernelU(cutoff, buffer1);
	v4CurrentValue.multAdd(-p(cutoff, cutoff + u), buffer1);
	if (u > 2.0f * cutoff)
	{
		integralKernelU(-cutoff, buffer1);
		v4CurrentValue.multAdd(-p(cutoff, cutoff - u), buffer1);
	}

	//Katanin contribution
	std::function<void(float, ValueSuperbundle<float, 2> &)> integralKernelSKatanin = [&](float wp, ValueSuperbundle<float, 2> &returnBuffer)->void { integralKernelS(wp, returnBuffer); returnBuffer *= pKataninContribution(wp, s + wp); };
	std::function<void(float, ValueSuperbundle<float, 2> &)> integralKernelTKatanin = [&](float wp, ValueSuperbundle<float, 2> &returnBuffer)->void { integralKernelT(wp, returnBuffer); returnBuffer *= pKataninContribution(wp, t + wp); };
	std::function<void(float, ValueSuperbundle<float, 2> &)> integralKernelUKatanin = [&](float wp, ValueSuperbundle<float, 2> &returnBuffer)->void { integralKernelU(wp, returnBuffer); returnBuffer *= -pKataninContribution(wp, u + wp); };

	if (-(s + cutoff) > *FrgCommon::frequency().beginNegative())
	{
		ImplicitIntegrator::integrateWithObscureRightBoundary(FrgCommon::frequency().beginNegative(), -(s + cutoff), integralKernelSKatanin, buffer1, buffer2);
		v4CurrentValue += buffer2;
	}
	if (s - cutoff > cutoff)
	{
		ImplicitIntegrator::integrateWithObscureBoundaries(cutoff - s, -cutoff, integralKernelSKatanin, buffer1, buffer2);
		v4CurrentValue += buffer2;
	}
	if (cutoff < *FrgCommon::frequency().last())
	{
		ImplicitIntegrator::integrateWithObscureLeftBoundary(cutoff, FrgCommon::frequency().last(), integralKernelSKatanin, buffer1, buffer2);
		v4CurrentValue += buffer2;
	}

	if (-(t + cutoff) > *FrgCommon::frequency().beginNegative())
	{
		ImplicitIntegrator::integrateWithObscureRightBoundary(FrgCommon::frequency().beginNegative(), -(t + cutoff), integralKernelTKatanin, buffer1, buffer2);
		v4CurrentValue += buffer2;
	}
	if (t - cutoff > cutoff)
	{
		ImplicitIntegrator::integrateWithObscureBoundaries(cutoff - t, -cutoff, integralKernelTKatanin, buffer1, buffer2);
		v4CurrentValue += buffer2;
	}
	if (cutoff < *FrgCommon::frequency().last())
	{
		ImplicitIntegrator::integrateWithObscureLeftBoundary(cutoff, FrgCommon::frequency().last(), integralKernelTKatanin, buffer1, buffer2);
		v4CurrentValue += buffer2;
	}

	if (-(u + cutoff) > *FrgCommon::frequency().beginNegative())
	{
		ImplicitIntegrator::integrateWithObscureRightBoundary(FrgCommon::frequency().beginNegative(), -(u + cutoff), integralKernelUKatanin, buffer1, buffer2);
		v4CurrentValue += buffer2;
	}
	if (u - cutoff > cutoff)
	{
		ImplicitIntegrator::integrateWithObscureBoundaries(cutoff - u, -cutoff, integralKernelUKatanin, buffer1, buffer2);
		v4CurrentValue += buffer2;
	}
	if (cutoff < *FrgCommon::frequency().last())
	{
		ImplicitIntegrator::integrateWithObscureLeftBoundary(cutoff, FrgCommon::frequency().last(), integralKernelUKatanin, buffer1, buffer2);
		v4CurrentValue += buffer2;
	}

	//prefactor
	v4CurrentValue /= 2.0f * (float)M_PI;

	for (int rid = 0; rid < FrgCommon::lattice().size; ++rid) static_cast<SU2EffectiveAction *>(_flow)->vertexTwoParticle->_dataSS[iterator * FrgCommon::lattice().size + rid] = v4CurrentValue.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Spin))[rid];
	for (int rid = 0; rid < FrgCommon::lattice().size; ++rid) static_cast<SU2EffectiveAction *>(_flow)->vertexTwoParticle->_dataDD[iterator * FrgCommon::lattice().size + rid] = v4CurrentValue.bundle(static_cast<int>(SU2VertexTwoParticle::Symmetry::Density))[rid];
}