/**
 * @file TRIFrgCore.cpp
 * @author Finn Lasse Buessen
 * @brief FrgCore implementation for time reversal invariant models.
 * 
 * @copyright Copyright (c) 2020
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include "lib/InputParser.hpp"
#include "lib/Integrator.hpp"
#include "SpinParser.hpp"
#include "TRIFrgCore.hpp"
#include "TRIEffectiveAction.hpp"

TRIFrgCore::TRIFrgCore(const SpinModel &spinModel, const std::vector<Measurement *> &measurements, const std::map<std::string, std::string> &options) : FrgCore(measurements)
{
	//init options
	normalization = NAN;

	for (auto option : options)
	{
		if (option.first == "normalization") normalization = InputParser::stringToFloat(option.second);
		else throw Exception(Exception::Type::InitializationError, "Unknown spin model option '" + option.first + "'.");
	}
	if (std::isnan(normalization)) normalization = 1.0f;

	Log::log << Log::LogLevel::Info << "FRG core energy normalization is set to " << normalization << "." << Log::endl;

	//init data
	_flowingFunctional = new TRIEffectiveAction(*FrgCommon::cutoff().begin(), spinModel, this);
	_flow = new TRIEffectiveAction();

	//init loadManager
	//stack0
	dataStacks[0] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		&_flowingFunctional->cutoff,
		1);
	//stack1
	dataStacks[1] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexSingleParticle->_data,
		static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexSingleParticle->size);
	//stack2
	dataStacks[2] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_data,
		static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size);
	//stack3
	dataStacks[3] = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		&_flow->cutoff,
		1,
		[&](int) { _flow->cutoff = static_cast<TRIEffectiveAction *>(_flowingFunctional)->cutoff; });
	//stack4
	dataStacks[4] = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		static_cast<TRIEffectiveAction *>(_flow)->vertexSingleParticle->_data,
		static_cast<TRIEffectiveAction *>(_flow)->vertexSingleParticle->size,
		[&](int x) { _calculateVertexSingleParticle(x); },
		1,
		1,
		1);
	//stack5
	dataStacks[5] = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		static_cast<TRIEffectiveAction *>(_flow)->vertexTwoParticle->_data,
		static_cast<TRIEffectiveAction *>(_flow)->vertexTwoParticle->sizeFrequency,
		[&](int x) { _calculateVertexTwoParticle(x); },
		16 * FrgCommon::lattice().size,
		FrgCommon::frequency().size);
}

TRIFrgCore::~TRIFrgCore()
{
	delete _flowingFunctional;
	delete _flow;
}

void TRIFrgCore::computeStep()
{
	//update cutoff and broadcast
	SpinParser::spinParser()->getLoadManager()->calculate(dataStacks[3]);
	SpinParser::spinParser()->getLoadManager()->broadcast(dataStacks[3]);
	//calculate 1-particle vertices and broadcast (required for Katanin calculation)
	SpinParser::spinParser()->getLoadManager()->calculate(dataStacks[4]);
	SpinParser::spinParser()->getLoadManager()->broadcast(dataStacks[4]);
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
	managedMeasurementStacks.push_back(dataStacks[5]);
	SpinParser::spinParser()->getLoadManager()->calculate(managedMeasurementStacks.data(), int(managedMeasurementStacks.size()));
}

void TRIFrgCore::finalizeStep(float newCutoff)
{
	//determine cutoff set
	float cutoffStep = newCutoff - _flowingFunctional->cutoff;

	//set new cutoff value
	_flowingFunctional->cutoff = newCutoff;

	//add _flow to single particle vertex
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexSingleParticle->size; ++i) static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexSingleParticle->_data[i] += cutoffStep * static_cast<TRIEffectiveAction *>(_flow)->vertexSingleParticle->_data[i];

	//add _flow to two particle vertex
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size; ++i) static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_data[i] += cutoffStep * static_cast<TRIEffectiveAction *>(_flow)->vertexTwoParticle->_data[i];

	//broadcast updated effective action
	SpinParser::spinParser()->getLoadManager()->broadcast({ dataStacks[0], dataStacks[1], dataStacks[2] });
}

void TRIFrgCore::_calculateVertexSingleParticle(const int iterator)
{
	float cutoff = _flowingFunctional->cutoff;
	TRIVertexSingleParticle *v2 = static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexSingleParticle;
	TRIVertexTwoParticle *v4 = static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexTwoParticle;

	float w;
	v2->expandIterator(iterator, w);
	float v2CurrentValue = 0;

	//term1
	float sum = 0;
	for (auto j = FrgCommon::lattice().getRange(0); j != FrgCommon::lattice().end(); ++j)
	{
		sum += v4->getValue(FrgCommon::lattice().zero(), j, w + cutoff, 0.0f, w - cutoff, SpinComponent::None, SpinComponent::None, TRIVertexTwoParticle::FrequencyChannel::None);
		sum -= v4->getValue(FrgCommon::lattice().zero(), j, w - cutoff, 0.0f, w + cutoff, SpinComponent::None, SpinComponent::None, TRIVertexTwoParticle::FrequencyChannel::None);
	}
	v2CurrentValue -= 2.0f * sum;

	//term2
	v2CurrentValue += v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SpinComponent::X, SpinComponent::X, TRIVertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SpinComponent::X, SpinComponent::X, TRIVertexTwoParticle::FrequencyChannel::None);
	v2CurrentValue += v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SpinComponent::Y, SpinComponent::Y, TRIVertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SpinComponent::Y, SpinComponent::Y, TRIVertexTwoParticle::FrequencyChannel::None);
	v2CurrentValue += v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SpinComponent::Z, SpinComponent::Z, TRIVertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SpinComponent::Z, SpinComponent::Z, TRIVertexTwoParticle::FrequencyChannel::None);
	v2CurrentValue += v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SpinComponent::None, SpinComponent::None, TRIVertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SpinComponent::None, SpinComponent::None, TRIVertexTwoParticle::FrequencyChannel::None);

	//prefactor
	v2CurrentValue /= (2.0f * (float)M_PI * (cutoff + v2->getValue(cutoff)));

	static_cast<TRIEffectiveAction *>(_flow)->vertexSingleParticle->_data[iterator] = v2CurrentValue;
}

void TRIFrgCore::_calculateVertexTwoParticle(const int iterator)
{
	float cutoff = _flowingFunctional->cutoff;
	TRIVertexSingleParticle *v2 = static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexSingleParticle;
	TRIVertexTwoParticle *v4 = static_cast<TRIEffectiveAction *>(_flowingFunctional)->vertexTwoParticle;

	float s, t, u;
	v4->expandIterator(iterator, s, t, u);

	//vertex buffers
	ValueSuperbundle<float, 16> buffer1(FrgCommon::lattice().size);
	ValueSuperbundle<float, 16> buffer2(FrgCommon::lattice().size);
	ValueSuperbundle<float, 16> bufferRPA(FrgCommon::lattice().size);
	ValueSuperbundle<float, 16> v4CurrentValue(FrgCommon::lattice().size);
	ValueSuperbundle<float, 16> stackBuffers[4] = {
		ValueSuperbundle<float, 16>(FrgCommon::lattice().size),
		ValueSuperbundle<float, 16>(FrgCommon::lattice().size),
		ValueSuperbundle<float, 16>(FrgCommon::lattice().size),
		ValueSuperbundle<float, 16>(FrgCommon::lattice().size)
	};

	//transfer frequencies
	float w1p = 0.5f * (s + t + u);
	float w1 = 0.5f * (s - t + u);
	float w2p = 0.5f * (s - t - u);
	float w2 = 0.5f * (s + t - u);

	//integrand of the ferquency integral
	auto integralKernelS = [&](const float wp, ValueSuperbundle<float, 16> &returnBuffer) -> void
	{
		//pp-ladder A and B (positive sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab0 = v4->generateAccessBuffer(s, w2 + wp, w1 + wp, TRIVertexTwoParticle::FrequencyChannel::S);
		const TRIVertexTwoParticleAccessBuffer<4> ab1 = v4->generateAccessBuffer(s, -w2p - wp, w1p + wp, TRIVertexTwoParticle::FrequencyChannel::S);
		//pp-ladder A and B (positive sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab2 = v4->generateAccessBuffer(s, -w1 - wp, -w2 - wp, TRIVertexTwoParticle::FrequencyChannel::S);
		const TRIVertexTwoParticleAccessBuffer<4> ab3 = v4->generateAccessBuffer(s, w1p + wp, -w2p - wp, TRIVertexTwoParticle::FrequencyChannel::S);

		v4->getValueSuperbundle(ab0, stackBuffers[0]);
		v4->getValueSuperbundle(ab1, stackBuffers[1]);
		v4->getValueSuperbundle(ab2, stackBuffers[2]);
		v4->getValueSuperbundle(ab3, stackBuffers[3]);

		//calculate _flow
		returnBuffer.reset();

		#pragma region ppLadder
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(15));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(15));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(14));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(14));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(13));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(13));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(12));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(12));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(11));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(11));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(10));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(10));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(9));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(9));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(8));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(8));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(7));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(7));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(6));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(6));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(5));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(5));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(4));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(4));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(3));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(3));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(2));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(2));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(1));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(1));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(0));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(0));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(12));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(12));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(13));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(13));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(14));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(14));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(15));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(15));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(8));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(8));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(9));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(9));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(10));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(10));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(11));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(11));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(4));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(4));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(5));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(5));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(6));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(6));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(7));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(7));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(0));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(0));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(1));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(1));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(2));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(2));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(3));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(3));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(13));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(13));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(12));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(12));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(15));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(15));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(14));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(14));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(9));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(9));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(8));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(8));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(11));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(11));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(10));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(10));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(5));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(5));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(4));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(4));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(7));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(7));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(6));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(6));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(1));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(1));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(0));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(0));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(3));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(3));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(2));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(2));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(14));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(14));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(15));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(15));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(12));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(12));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(13));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(13));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(10));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(10));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(11));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(11));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(8));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(8));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(9));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(9));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(6));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(6));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(7));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(7));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(4));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(4));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(5));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(5));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(2));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(2));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(3));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(3));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(0));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(0));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(1));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(1));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(3));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(3));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(2));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(2));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(1));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(1));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(0));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(0));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(7));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(7));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(6));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(6));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(5));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(5));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(4));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(4));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(11));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(11));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(10));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(10));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(9));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(9));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(8));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(8));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(15));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(15));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(14));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(14));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(13));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(13));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(12));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(12));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(0));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(0));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(1));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(1));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(2));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(2));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(3));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(3));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(4));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(4));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(5));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(5));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(6));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(6));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(7));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(7));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(8));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(8));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(9));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(9));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(10));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(10));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(11));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(11));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(12));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(12));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(13));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(13));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(14));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(14));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(15));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(15));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(1));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(1));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(0));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(0));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(3));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(3));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(2));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(2));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(5));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(5));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(4));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(4));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(7));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(7));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(6));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(6));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(9));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(9));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(8));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(8));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(11));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(11));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(10));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(10));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(13));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(13));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(12));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(12));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(15));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(15));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(14));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(14));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(2));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(2));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(3));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(3));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(0));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(0));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(1));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(1));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(6));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(6));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(7));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(7));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(4));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(4));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(5));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(5));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(10));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(10));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(11));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(11));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(8));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(8));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(9));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(9));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(14));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(14));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(15));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(15));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(12));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(12));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(13));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(13));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(7));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(7));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(6));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(6));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(5));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(5));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(4));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(4));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(3));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(3));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(2));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(2));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(1));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(1));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(0));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(0));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(15));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(15));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(14));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(14));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(13));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(13));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(12));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(12));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(11));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(11));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(10));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(10));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(9));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(9));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(8));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(8));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(4));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(4));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(5));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(5));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(6));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(6));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(7));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(7));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(0));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(0));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(1));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(1));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(2));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(2));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(3));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(3));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(12));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(12));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(13));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(13));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(14));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(14));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(15));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(15));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(8));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(8));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(9));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(9));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(10));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(10));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(11));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(11));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(5));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(5));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(4));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(4));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(7));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(7));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(6));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(6));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(1));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(1));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(0));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(0));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(3));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(3));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(2));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(2));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(13));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(13));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(12));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(12));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(15));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(15));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(14));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(14));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(9));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(9));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(8));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(8));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(11));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(11));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(10));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(10));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(6));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(6));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(7));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(7));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(4));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(4));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(5));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(5));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(2));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(2));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(3));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(3));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(0));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(0));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(1));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(1));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(14));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(14));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(15));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(15));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(12));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(12));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(13));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(13));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(10));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(10));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(11));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(11));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(8));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(8));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(9));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(9));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(11));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(11));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(10));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(10));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(9));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(9));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(8));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(8));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(15));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(15));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(14));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(14));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(13));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(13));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(12));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(12));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(3));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(3));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(2));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(2));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(1));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(1));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(0));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(0));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(7));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(7));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(6));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(6));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(5));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(5));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(4));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(4));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(8));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(8));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(9));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(9));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(10));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(10));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(11));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(11));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(12));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(12));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(13));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(13));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(14));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(14));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(15));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(15));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(0));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(0));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(1));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(1));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(2));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(2));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(3));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(3));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(4));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(4));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(5));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(5));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(6));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(6));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(7));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(7));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(9));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(9));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(8));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(8));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(11));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(11));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(10));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(10));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(13));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(13));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(12));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(12));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(15));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(15));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(14));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(14));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(1));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(1));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(0));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(0));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(3));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(3));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(2));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(2));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(5));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(5));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(4));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(4));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(7));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(7));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(6));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(6));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(15), stackBuffers[1].bundle(10));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(15), stackBuffers[3].bundle(10));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(11));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(11));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(8));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(8));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(9));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(9));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(14));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(14));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(15));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(15));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(12));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(12));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(13));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(13));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(2));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(2));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(3));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(3));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(0));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(0));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(1));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(1));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(6));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(6));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(7));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(7));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(4));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(4));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(5));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(5));
		#pragma endregion
	};

	auto integralKernelT = [&](const float wp, ValueSuperbundle<float, 16> &returnBuffer) -> void
	{
		//RPA diagram A and B equal chalice diagram A and inverse chalice diagram B, respectively (negative sign)
		//chalice diagram A (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab0 = v4->generateAccessBuffer(w1 - wp, t, w1p + wp, TRIVertexTwoParticle::FrequencyChannel::T);
		//inverse chalice diagram B (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab1 = v4->generateAccessBuffer(w2p - wp, t, -w2 - wp, TRIVertexTwoParticle::FrequencyChannel::T);
		//chalice diagram A (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab2 = v4->generateAccessBuffer(w1p + wp, t, w1 - wp, TRIVertexTwoParticle::FrequencyChannel::T);
		//inverse chalice diagram B (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab3 = v4->generateAccessBuffer(w2 + wp, t, -w2p + wp, TRIVertexTwoParticle::FrequencyChannel::T);

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

			#pragma region RPA
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(15)[rid] += 2 * stackBuffers[0].bundle(15).data()[overlap.rid1[i]] * stackBuffers[1].bundle(15).data()[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(15)[rid] += 2 * stackBuffers[2].bundle(15).data()[overlap.rid1[i]] * stackBuffers[3].bundle(15).data()[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(15)[rid] -= 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(15)[rid] -= 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(15)[rid] -= 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(15)[rid] -= 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(15)[rid] -= 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(15)[rid] -= 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(12)[rid] += 2 * stackBuffers[0].bundle(15).data()[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(12)[rid] += 2 * stackBuffers[2].bundle(15).data()[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(12)[rid] += 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(12)[rid] += 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(12)[rid] += 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(12)[rid] += 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(12)[rid] += 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(12)[rid] += 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(13)[rid] += 2 * stackBuffers[0].bundle(15).data()[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(13)[rid] += 2 * stackBuffers[2].bundle(15).data()[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(13)[rid] += 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(13)[rid] += 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(13)[rid] += 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(13)[rid] += 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(13)[rid] += 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(13)[rid] += 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(14)[rid] += 2 * stackBuffers[0].bundle(15).data()[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(14)[rid] += 2 * stackBuffers[2].bundle(15).data()[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(14)[rid] += 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(14)[rid] += 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(14)[rid] += 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(14)[rid] += 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(14)[rid] += 2 * stackBuffers[0].bundle(12 + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(14)[rid] += 2 * stackBuffers[2].bundle(12 + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(3)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(15).data()[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(3)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(15).data()[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(3)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(3)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(3)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(3)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(3)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(3)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(0)[rid] -= 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(0)[rid] -= 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(0)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(0)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(0)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(0)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(0)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(0)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(1)[rid] -= 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(1)[rid] -= 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(1)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(1)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(1)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(1)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(1)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(1)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(2)[rid] -= 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(2)[rid] -= 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(2)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(2)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(2)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(2)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(2)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(2)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedX1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(7)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(15).data()[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(7)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(15).data()[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(7)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(7)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(7)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(7)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(7)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(7)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(4)[rid] -= 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(4)[rid] -= 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(4)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(4)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(4)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(4)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(4)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(4)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(5)[rid] -= 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(5)[rid] -= 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(5)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(5)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(5)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(5)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(5)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(5)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(6)[rid] -= 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(6)[rid] -= 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(6)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(6)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(6)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(6)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(6)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(6)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedY1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(11)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(15).data()[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(11)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(15).data()[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(11)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(11)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(11)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(11)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(11)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(11)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + 3)[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(8)[rid] -= 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(8)[rid] -= 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(8)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(8)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(8)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(8)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(8)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(8)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(9)[rid] -= 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(9)[rid] -= 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(9)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(9)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(9)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(9)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(9)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(9)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(10)[rid] -= 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + 3)[overlap.rid1[i]] * stackBuffers[1].bundle(12 + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(10)[rid] -= 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + 3)[overlap.rid1[i]] * stackBuffers[3].bundle(12 + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(10)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(10)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedZ2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(10)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(10)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedY2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(10)[rid] += 2 * stackBuffers[0].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(10)[rid] += 2 * stackBuffers[2].bundle(4 * static_cast<int>(overlap.transformedZ1[i]) + static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[3].bundle(4 * static_cast<int>(overlap.transformedX2[i]) + static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}
			#pragma endregion
		}
		returnBuffer += bufferRPA;

		//chalice diagram B (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab4 = v4->generateAccessBuffer(w2p - wp, -w2 - wp, t, TRIVertexTwoParticle::FrequencyChannel::U);
		//chalice diagram B (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab5 = v4->generateAccessBuffer(w2 + wp, -w2p + wp, t, TRIVertexTwoParticle::FrequencyChannel::U);

		const float valLocal4[16] = {
			v4->getValueLocal(SpinComponent::X, SpinComponent::X, ab4),
			v4->getValueLocal(SpinComponent::X, SpinComponent::Y, ab4),
			v4->getValueLocal(SpinComponent::X, SpinComponent::Z, ab4),
			v4->getValueLocal(SpinComponent::X, SpinComponent::None, ab4),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::X, ab4),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::Y, ab4),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::Z, ab4),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::None, ab4),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::X, ab4),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::Y, ab4),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::Z, ab4),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::None, ab4),
			v4->getValueLocal(SpinComponent::None, SpinComponent::X, ab4),
			v4->getValueLocal(SpinComponent::None, SpinComponent::Y, ab4),
			v4->getValueLocal(SpinComponent::None, SpinComponent::Z, ab4),
			v4->getValueLocal(SpinComponent::None, SpinComponent::None, ab4)
		};
		const float valLocal5[16] = {
			v4->getValueLocal(SpinComponent::X, SpinComponent::X, ab5),
			v4->getValueLocal(SpinComponent::X, SpinComponent::Y, ab5),
			v4->getValueLocal(SpinComponent::X, SpinComponent::Z, ab5),
			v4->getValueLocal(SpinComponent::X, SpinComponent::None, ab5),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::X, ab5),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::Y, ab5),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::Z, ab5),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::None, ab5),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::X, ab5),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::Y, ab5),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::Z, ab5),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::None, ab5),
			v4->getValueLocal(SpinComponent::None, SpinComponent::X, ab5),
			v4->getValueLocal(SpinComponent::None, SpinComponent::Y, ab5),
			v4->getValueLocal(SpinComponent::None, SpinComponent::Z, ab5),
			v4->getValueLocal(SpinComponent::None, SpinComponent::None, ab5)
		};

		#pragma region chalice
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(15), valLocal4[15]);
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(15), valLocal5[15]);
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(15), valLocal4[10]);
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(15), valLocal5[10]);
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(15), valLocal4[5]);
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(15), valLocal5[5]);
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(15), valLocal4[0]);
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(15), valLocal5[0]);
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(14), valLocal4[14]);
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(14), valLocal5[14]);
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(14), valLocal4[11]);
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(14), valLocal5[11]);
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(14), valLocal4[4]);
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(14), valLocal5[4]);
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(14), valLocal4[1]);
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(14), valLocal5[1]);
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(13), valLocal4[13]);
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(13), valLocal5[13]);
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(13), valLocal4[8]);
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(13), valLocal5[8]);
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(13), valLocal4[7]);
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(13), valLocal5[7]);
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(13), valLocal4[2]);
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(13), valLocal5[2]);
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(12), valLocal4[12]);
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(12), valLocal5[12]);
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(12), valLocal4[9]);
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(12), valLocal5[9]);
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(12), valLocal4[6]);
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(12), valLocal5[6]);
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(12), valLocal4[3]);
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(12), valLocal5[3]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(15), valLocal4[12]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(15), valLocal5[12]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(15), valLocal4[9]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(15), valLocal5[9]);
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(15), valLocal4[6]);
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(15), valLocal5[6]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(15), valLocal4[3]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(15), valLocal5[3]);
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(14), valLocal4[13]);
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(14), valLocal5[13]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(14), valLocal4[8]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(14), valLocal5[8]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(14), valLocal4[7]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(14), valLocal5[7]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(14), valLocal4[2]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(14), valLocal5[2]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(13), valLocal4[14]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(13), valLocal5[14]);
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(13), valLocal4[11]);
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(13), valLocal5[11]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(13), valLocal4[4]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(13), valLocal5[4]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(13), valLocal4[1]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(13), valLocal5[1]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(12), valLocal4[15]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(12), valLocal5[15]);
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(12), valLocal4[10]);
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(12), valLocal5[10]);
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(12), valLocal4[5]);
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(12), valLocal5[5]);
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(12), valLocal4[0]);
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(12), valLocal5[0]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(15), valLocal4[13]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(15), valLocal5[13]);
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(15), valLocal4[8]);
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(15), valLocal5[8]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(15), valLocal4[7]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(15), valLocal5[7]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(15), valLocal4[2]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(15), valLocal5[2]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(14), valLocal4[12]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(14), valLocal5[12]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(14), valLocal4[9]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(14), valLocal5[9]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(14), valLocal4[6]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(14), valLocal5[6]);
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(14), valLocal4[3]);
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(14), valLocal5[3]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(13), valLocal4[15]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(13), valLocal5[15]);
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(13), valLocal4[10]);
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(13), valLocal5[10]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(13), valLocal4[5]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(13), valLocal5[5]);
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(13), valLocal4[0]);
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(13), valLocal5[0]);
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(12), valLocal4[14]);
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(12), valLocal5[14]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(12), valLocal4[11]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(12), valLocal5[11]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(12), valLocal4[4]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(12), valLocal5[4]);
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(12), valLocal4[1]);
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(12), valLocal5[1]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(15), valLocal4[14]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(15), valLocal5[14]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(15), valLocal4[11]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(15), valLocal5[11]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(15), valLocal4[4]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(15), valLocal5[4]);
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(15), valLocal4[1]);
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(15), valLocal5[1]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(14), valLocal4[15]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(14), valLocal5[15]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(14), valLocal4[10]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(14), valLocal5[10]);
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(14), valLocal4[5]);
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(14), valLocal5[5]);
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(14), valLocal4[0]);
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(14), valLocal5[0]);
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(13), valLocal4[12]);
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(13), valLocal5[12]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(13), valLocal4[9]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(13), valLocal5[9]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(13), valLocal4[6]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(13), valLocal5[6]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(13), valLocal4[3]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(13), valLocal5[3]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(12), valLocal4[13]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(12), valLocal5[13]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(12), valLocal4[8]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(12), valLocal5[8]);
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(12), valLocal4[7]);
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(12), valLocal5[7]);
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(12), valLocal4[2]);
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(12), valLocal5[2]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(3), valLocal4[15]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(3), valLocal5[15]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(3), valLocal4[10]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(3), valLocal5[10]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(3), valLocal4[5]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(3), valLocal5[5]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(3), valLocal4[0]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(3), valLocal5[0]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(2), valLocal4[14]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(2), valLocal5[14]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(2), valLocal4[11]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(2), valLocal5[11]);
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(2), valLocal4[4]);
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(2), valLocal5[4]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(2), valLocal4[1]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(2), valLocal5[1]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(1), valLocal4[13]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(1), valLocal5[13]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(1), valLocal4[8]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(1), valLocal5[8]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(1), valLocal4[7]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(1), valLocal5[7]);
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(1), valLocal4[2]);
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(1), valLocal5[2]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(0), valLocal4[12]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(0), valLocal5[12]);
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(0), valLocal4[9]);
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(0), valLocal5[9]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(0), valLocal4[6]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(0), valLocal5[6]);
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(0), valLocal4[3]);
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(0), valLocal5[3]);
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(3), valLocal4[12]);
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(3), valLocal5[12]);
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(3), valLocal4[9]);
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(3), valLocal5[9]);
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(3), valLocal4[6]);
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(3), valLocal5[6]);
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(3), valLocal4[3]);
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(3), valLocal5[3]);
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(2), valLocal4[13]);
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(2), valLocal5[13]);
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(2), valLocal4[8]);
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(2), valLocal5[8]);
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(2), valLocal4[7]);
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(2), valLocal5[7]);
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(2), valLocal4[2]);
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(2), valLocal5[2]);
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(1), valLocal4[14]);
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(1), valLocal5[14]);
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(1), valLocal4[11]);
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(1), valLocal5[11]);
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(1), valLocal4[4]);
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(1), valLocal5[4]);
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(1), valLocal4[1]);
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(1), valLocal5[1]);
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(0), valLocal4[15]);
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(0), valLocal5[15]);
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(0), valLocal4[10]);
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(0), valLocal5[10]);
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(0), valLocal4[5]);
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(0), valLocal5[5]);
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(0), valLocal4[0]);
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(0), valLocal5[0]);
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(3), valLocal4[13]);
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(3), valLocal5[13]);
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(3), valLocal4[8]);
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(3), valLocal5[8]);
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(3), valLocal4[7]);
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(3), valLocal5[7]);
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(3), valLocal4[2]);
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(3), valLocal5[2]);
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(2), valLocal4[12]);
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(2), valLocal5[12]);
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(2), valLocal4[9]);
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(2), valLocal5[9]);
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(2), valLocal4[6]);
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(2), valLocal5[6]);
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(2), valLocal4[3]);
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(2), valLocal5[3]);
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(1), valLocal4[15]);
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(1), valLocal5[15]);
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(1), valLocal4[10]);
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(1), valLocal5[10]);
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(1), valLocal4[5]);
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(1), valLocal5[5]);
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(1), valLocal4[0]);
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(1), valLocal5[0]);
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(0), valLocal4[14]);
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(0), valLocal5[14]);
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(0), valLocal4[11]);
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(0), valLocal5[11]);
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(0), valLocal4[4]);
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(0), valLocal5[4]);
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(0), valLocal4[1]);
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(0), valLocal5[1]);
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(3), valLocal4[14]);
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(3), valLocal5[14]);
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(3), valLocal4[11]);
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(3), valLocal5[11]);
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(3), valLocal4[4]);
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(3), valLocal5[4]);
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(3), valLocal4[1]);
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(3), valLocal5[1]);
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(2), valLocal4[15]);
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(2), valLocal5[15]);
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(2), valLocal4[10]);
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(2), valLocal5[10]);
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(2), valLocal4[5]);
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(2), valLocal5[5]);
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(2), valLocal4[0]);
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(2), valLocal5[0]);
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(1), valLocal4[12]);
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(1), valLocal5[12]);
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(1), valLocal4[9]);
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(1), valLocal5[9]);
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(1), valLocal4[6]);
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(1), valLocal5[6]);
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(1), valLocal4[3]);
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(1), valLocal5[3]);
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(0), valLocal4[13]);
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(0), valLocal5[13]);
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(0), valLocal4[8]);
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(0), valLocal5[8]);
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(0), valLocal4[7]);
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(0), valLocal5[7]);
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(0), valLocal4[2]);
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(0), valLocal5[2]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(7), valLocal4[15]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(7), valLocal5[15]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(7), valLocal4[10]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(7), valLocal5[10]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(7), valLocal4[5]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(7), valLocal5[5]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(7), valLocal4[0]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(7), valLocal5[0]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(6), valLocal4[14]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(6), valLocal5[14]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(6), valLocal4[11]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(6), valLocal5[11]);
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(6), valLocal4[4]);
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(6), valLocal5[4]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(6), valLocal4[1]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(6), valLocal5[1]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(5), valLocal4[13]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(5), valLocal5[13]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(5), valLocal4[8]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(5), valLocal5[8]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(5), valLocal4[7]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(5), valLocal5[7]);
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(5), valLocal4[2]);
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(5), valLocal5[2]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(4), valLocal4[12]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(4), valLocal5[12]);
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(4), valLocal4[9]);
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(4), valLocal5[9]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(4), valLocal4[6]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(4), valLocal5[6]);
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(4), valLocal4[3]);
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(4), valLocal5[3]);
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(7), valLocal4[12]);
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(7), valLocal5[12]);
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(7), valLocal4[9]);
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(7), valLocal5[9]);
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(7), valLocal4[6]);
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(7), valLocal5[6]);
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(7), valLocal4[3]);
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(7), valLocal5[3]);
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(6), valLocal4[13]);
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(6), valLocal5[13]);
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(6), valLocal4[8]);
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(6), valLocal5[8]);
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(6), valLocal4[7]);
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(6), valLocal5[7]);
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(6), valLocal4[2]);
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(6), valLocal5[2]);
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(5), valLocal4[14]);
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(5), valLocal5[14]);
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(5), valLocal4[11]);
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(5), valLocal5[11]);
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(5), valLocal4[4]);
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(5), valLocal5[4]);
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(5), valLocal4[1]);
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(5), valLocal5[1]);
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(4), valLocal4[15]);
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(4), valLocal5[15]);
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(4), valLocal4[10]);
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(4), valLocal5[10]);
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(4), valLocal4[5]);
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(4), valLocal5[5]);
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(4), valLocal4[0]);
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(4), valLocal5[0]);
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(7), valLocal4[13]);
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(7), valLocal5[13]);
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(7), valLocal4[8]);
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(7), valLocal5[8]);
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(7), valLocal4[7]);
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(7), valLocal5[7]);
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(7), valLocal4[2]);
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(7), valLocal5[2]);
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(6), valLocal4[12]);
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(6), valLocal5[12]);
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(6), valLocal4[9]);
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(6), valLocal5[9]);
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(6), valLocal4[6]);
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(6), valLocal5[6]);
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(6), valLocal4[3]);
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(6), valLocal5[3]);
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(5), valLocal4[15]);
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(5), valLocal5[15]);
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(5), valLocal4[10]);
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(5), valLocal5[10]);
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(5), valLocal4[5]);
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(5), valLocal5[5]);
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(5), valLocal4[0]);
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(5), valLocal5[0]);
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(4), valLocal4[14]);
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(4), valLocal5[14]);
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(4), valLocal4[11]);
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(4), valLocal5[11]);
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(4), valLocal4[4]);
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(4), valLocal5[4]);
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(4), valLocal4[1]);
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(4), valLocal5[1]);
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(7), valLocal4[14]);
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(7), valLocal5[14]);
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(7), valLocal4[11]);
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(7), valLocal5[11]);
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(7), valLocal4[4]);
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(7), valLocal5[4]);
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(7), valLocal4[1]);
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(7), valLocal5[1]);
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(6), valLocal4[15]);
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(6), valLocal5[15]);
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(6), valLocal4[10]);
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(6), valLocal5[10]);
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(6), valLocal4[5]);
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(6), valLocal5[5]);
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(6), valLocal4[0]);
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(6), valLocal5[0]);
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(5), valLocal4[12]);
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(5), valLocal5[12]);
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(5), valLocal4[9]);
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(5), valLocal5[9]);
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(5), valLocal4[6]);
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(5), valLocal5[6]);
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(5), valLocal4[3]);
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(5), valLocal5[3]);
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(4), valLocal4[13]);
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(4), valLocal5[13]);
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(4), valLocal4[8]);
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(4), valLocal5[8]);
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(4), valLocal4[7]);
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(4), valLocal5[7]);
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(4), valLocal4[2]);
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(4), valLocal5[2]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(11), valLocal4[15]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(11), valLocal5[15]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(11), valLocal4[10]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(11), valLocal5[10]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(11), valLocal4[5]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(11), valLocal5[5]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(11), valLocal4[0]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(11), valLocal5[0]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(10), valLocal4[14]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(10), valLocal5[14]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(10), valLocal4[11]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(10), valLocal5[11]);
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(10), valLocal4[4]);
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(10), valLocal5[4]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(10), valLocal4[1]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(10), valLocal5[1]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(9), valLocal4[13]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(9), valLocal5[13]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(9), valLocal4[8]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(9), valLocal5[8]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(9), valLocal4[7]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(9), valLocal5[7]);
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(9), valLocal4[2]);
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(9), valLocal5[2]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(8), valLocal4[12]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(8), valLocal5[12]);
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(8), valLocal4[9]);
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(8), valLocal5[9]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(8), valLocal4[6]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(8), valLocal5[6]);
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(8), valLocal4[3]);
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(8), valLocal5[3]);
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(11), valLocal4[12]);
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(11), valLocal5[12]);
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(11), valLocal4[9]);
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(11), valLocal5[9]);
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(11), valLocal4[6]);
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(11), valLocal5[6]);
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(11), valLocal4[3]);
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(11), valLocal5[3]);
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(10), valLocal4[13]);
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(10), valLocal5[13]);
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(10), valLocal4[8]);
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(10), valLocal5[8]);
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(10), valLocal4[7]);
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(10), valLocal5[7]);
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(10), valLocal4[2]);
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(10), valLocal5[2]);
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(9), valLocal4[14]);
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(9), valLocal5[14]);
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(9), valLocal4[11]);
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(9), valLocal5[11]);
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(9), valLocal4[4]);
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(9), valLocal5[4]);
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(9), valLocal4[1]);
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(9), valLocal5[1]);
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(8), valLocal4[15]);
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(8), valLocal5[15]);
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(8), valLocal4[10]);
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(8), valLocal5[10]);
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(8), valLocal4[5]);
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(8), valLocal5[5]);
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(8), valLocal4[0]);
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(8), valLocal5[0]);
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(11), valLocal4[13]);
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(11), valLocal5[13]);
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(11), valLocal4[8]);
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(11), valLocal5[8]);
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(11), valLocal4[7]);
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(11), valLocal5[7]);
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(11), valLocal4[2]);
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(11), valLocal5[2]);
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(10), valLocal4[12]);
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(10), valLocal5[12]);
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(10), valLocal4[9]);
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(10), valLocal5[9]);
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(10), valLocal4[6]);
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(10), valLocal5[6]);
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(10), valLocal4[3]);
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(10), valLocal5[3]);
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(9), valLocal4[15]);
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(9), valLocal5[15]);
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(9), valLocal4[10]);
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(9), valLocal5[10]);
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(9), valLocal4[5]);
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(9), valLocal5[5]);
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(9), valLocal4[0]);
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(9), valLocal5[0]);
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(8), valLocal4[14]);
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(8), valLocal5[14]);
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(8), valLocal4[11]);
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(8), valLocal5[11]);
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(8), valLocal4[4]);
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(8), valLocal5[4]);
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(8), valLocal4[1]);
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(8), valLocal5[1]);
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(11), valLocal4[14]);
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(11), valLocal5[14]);
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(11), valLocal4[11]);
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(11), valLocal5[11]);
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(11), valLocal4[4]);
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(11), valLocal5[4]);
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(11), valLocal4[1]);
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(11), valLocal5[1]);
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(10), valLocal4[15]);
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(10), valLocal5[15]);
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(10), valLocal4[10]);
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(10), valLocal5[10]);
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(10), valLocal4[5]);
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(10), valLocal5[5]);
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(10), valLocal4[0]);
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(10), valLocal5[0]);
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(9), valLocal4[12]);
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(9), valLocal5[12]);
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(9), valLocal4[9]);
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(9), valLocal5[9]);
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(9), valLocal4[6]);
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(9), valLocal5[6]);
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(9), valLocal4[3]);
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(9), valLocal5[3]);
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(8), valLocal4[13]);
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(8), valLocal5[13]);
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(8), valLocal4[8]);
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(8), valLocal5[8]);
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(8), valLocal4[7]);
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(8), valLocal5[7]);
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(8), valLocal4[2]);
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(8), valLocal5[2]);
		#pragma endregion

		//inverse chalice diagram A (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab6 = v4->generateAccessBuffer(w1 - wp, -w1p - wp, -t, TRIVertexTwoParticle::FrequencyChannel::U);
		//inverse chalice diagram A (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab7 = v4->generateAccessBuffer(w1p + wp, -w1 + wp, -t, TRIVertexTwoParticle::FrequencyChannel::U);

		const float valLocal6[16] = {
			v4->getValueLocal(SpinComponent::X, SpinComponent::X, ab6),
			v4->getValueLocal(SpinComponent::X, SpinComponent::Y, ab6),
			v4->getValueLocal(SpinComponent::X, SpinComponent::Z, ab6),
			v4->getValueLocal(SpinComponent::X, SpinComponent::None, ab6),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::X, ab6),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::Y, ab6),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::Z, ab6),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::None, ab6),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::X, ab6),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::Y, ab6),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::Z, ab6),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::None, ab6),
			v4->getValueLocal(SpinComponent::None, SpinComponent::X, ab6),
			v4->getValueLocal(SpinComponent::None, SpinComponent::Y, ab6),
			v4->getValueLocal(SpinComponent::None, SpinComponent::Z, ab6),
			v4->getValueLocal(SpinComponent::None, SpinComponent::None, ab6)
		};
		const float valLocal7[16] = {
			v4->getValueLocal(SpinComponent::X, SpinComponent::X, ab7),
			v4->getValueLocal(SpinComponent::X, SpinComponent::Y, ab7),
			v4->getValueLocal(SpinComponent::X, SpinComponent::Z, ab7),
			v4->getValueLocal(SpinComponent::X, SpinComponent::None, ab7),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::X, ab7),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::Y, ab7),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::Z, ab7),
			v4->getValueLocal(SpinComponent::Y, SpinComponent::None, ab7),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::X, ab7),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::Y, ab7),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::Z, ab7),
			v4->getValueLocal(SpinComponent::Z, SpinComponent::None, ab7),
			v4->getValueLocal(SpinComponent::None, SpinComponent::X, ab7),
			v4->getValueLocal(SpinComponent::None, SpinComponent::Y, ab7),
			v4->getValueLocal(SpinComponent::None, SpinComponent::Z, ab7),
			v4->getValueLocal(SpinComponent::None, SpinComponent::None, ab7)
		};

		#pragma region inverseChalice
		returnBuffer.bundle(15).multSub(valLocal6[15], stackBuffers[1].bundle(15));
		returnBuffer.bundle(15).multSub(valLocal7[15], stackBuffers[3].bundle(15));
		returnBuffer.bundle(15).multAdd(valLocal6[14], stackBuffers[1].bundle(11));
		returnBuffer.bundle(15).multAdd(valLocal7[14], stackBuffers[3].bundle(11));
		returnBuffer.bundle(15).multAdd(valLocal6[13], stackBuffers[1].bundle(7));
		returnBuffer.bundle(15).multAdd(valLocal7[13], stackBuffers[3].bundle(7));
		returnBuffer.bundle(15).multAdd(valLocal6[12], stackBuffers[1].bundle(3));
		returnBuffer.bundle(15).multAdd(valLocal7[12], stackBuffers[3].bundle(3));
		returnBuffer.bundle(15).multAdd(valLocal6[11], stackBuffers[1].bundle(11));
		returnBuffer.bundle(15).multAdd(valLocal7[11], stackBuffers[3].bundle(11));
		returnBuffer.bundle(15).multSub(valLocal6[10], stackBuffers[1].bundle(15));
		returnBuffer.bundle(15).multSub(valLocal7[10], stackBuffers[3].bundle(15));
		returnBuffer.bundle(15).multSub(valLocal6[9], stackBuffers[1].bundle(3));
		returnBuffer.bundle(15).multSub(valLocal7[9], stackBuffers[3].bundle(3));
		returnBuffer.bundle(15).multAdd(valLocal6[8], stackBuffers[1].bundle(7));
		returnBuffer.bundle(15).multAdd(valLocal7[8], stackBuffers[3].bundle(7));
		returnBuffer.bundle(15).multAdd(valLocal6[7], stackBuffers[1].bundle(7));
		returnBuffer.bundle(15).multAdd(valLocal7[7], stackBuffers[3].bundle(7));
		returnBuffer.bundle(15).multAdd(valLocal6[6], stackBuffers[1].bundle(3));
		returnBuffer.bundle(15).multAdd(valLocal7[6], stackBuffers[3].bundle(3));
		returnBuffer.bundle(15).multSub(valLocal6[5], stackBuffers[1].bundle(15));
		returnBuffer.bundle(15).multSub(valLocal7[5], stackBuffers[3].bundle(15));
		returnBuffer.bundle(15).multSub(valLocal6[4], stackBuffers[1].bundle(11));
		returnBuffer.bundle(15).multSub(valLocal7[4], stackBuffers[3].bundle(11));
		returnBuffer.bundle(15).multAdd(valLocal6[3], stackBuffers[1].bundle(3));
		returnBuffer.bundle(15).multAdd(valLocal7[3], stackBuffers[3].bundle(3));
		returnBuffer.bundle(15).multSub(valLocal6[2], stackBuffers[1].bundle(7));
		returnBuffer.bundle(15).multSub(valLocal7[2], stackBuffers[3].bundle(7));
		returnBuffer.bundle(15).multAdd(valLocal6[1], stackBuffers[1].bundle(11));
		returnBuffer.bundle(15).multAdd(valLocal7[1], stackBuffers[3].bundle(11));
		returnBuffer.bundle(15).multSub(valLocal6[0], stackBuffers[1].bundle(15));
		returnBuffer.bundle(15).multSub(valLocal7[0], stackBuffers[3].bundle(15));
		returnBuffer.bundle(12).multSub(valLocal6[15], stackBuffers[1].bundle(12));
		returnBuffer.bundle(12).multSub(valLocal7[15], stackBuffers[3].bundle(12));
		returnBuffer.bundle(12).multSub(valLocal6[14], stackBuffers[1].bundle(8));
		returnBuffer.bundle(12).multSub(valLocal7[14], stackBuffers[3].bundle(8));
		returnBuffer.bundle(12).multSub(valLocal6[13], stackBuffers[1].bundle(4));
		returnBuffer.bundle(12).multSub(valLocal7[13], stackBuffers[3].bundle(4));
		returnBuffer.bundle(12).multSub(valLocal6[12], stackBuffers[1].bundle(0));
		returnBuffer.bundle(12).multSub(valLocal7[12], stackBuffers[3].bundle(0));
		returnBuffer.bundle(12).multSub(valLocal6[11], stackBuffers[1].bundle(8));
		returnBuffer.bundle(12).multSub(valLocal7[11], stackBuffers[3].bundle(8));
		returnBuffer.bundle(12).multSub(valLocal6[10], stackBuffers[1].bundle(12));
		returnBuffer.bundle(12).multSub(valLocal7[10], stackBuffers[3].bundle(12));
		returnBuffer.bundle(12).multAdd(valLocal6[9], stackBuffers[1].bundle(0));
		returnBuffer.bundle(12).multAdd(valLocal7[9], stackBuffers[3].bundle(0));
		returnBuffer.bundle(12).multSub(valLocal6[8], stackBuffers[1].bundle(4));
		returnBuffer.bundle(12).multSub(valLocal7[8], stackBuffers[3].bundle(4));
		returnBuffer.bundle(12).multSub(valLocal6[7], stackBuffers[1].bundle(4));
		returnBuffer.bundle(12).multSub(valLocal7[7], stackBuffers[3].bundle(4));
		returnBuffer.bundle(12).multSub(valLocal6[6], stackBuffers[1].bundle(0));
		returnBuffer.bundle(12).multSub(valLocal7[6], stackBuffers[3].bundle(0));
		returnBuffer.bundle(12).multSub(valLocal6[5], stackBuffers[1].bundle(12));
		returnBuffer.bundle(12).multSub(valLocal7[5], stackBuffers[3].bundle(12));
		returnBuffer.bundle(12).multAdd(valLocal6[4], stackBuffers[1].bundle(8));
		returnBuffer.bundle(12).multAdd(valLocal7[4], stackBuffers[3].bundle(8));
		returnBuffer.bundle(12).multSub(valLocal6[3], stackBuffers[1].bundle(0));
		returnBuffer.bundle(12).multSub(valLocal7[3], stackBuffers[3].bundle(0));
		returnBuffer.bundle(12).multAdd(valLocal6[2], stackBuffers[1].bundle(4));
		returnBuffer.bundle(12).multAdd(valLocal7[2], stackBuffers[3].bundle(4));
		returnBuffer.bundle(12).multSub(valLocal6[1], stackBuffers[1].bundle(8));
		returnBuffer.bundle(12).multSub(valLocal7[1], stackBuffers[3].bundle(8));
		returnBuffer.bundle(12).multSub(valLocal6[0], stackBuffers[1].bundle(12));
		returnBuffer.bundle(12).multSub(valLocal7[0], stackBuffers[3].bundle(12));
		returnBuffer.bundle(13).multSub(valLocal6[15], stackBuffers[1].bundle(13));
		returnBuffer.bundle(13).multSub(valLocal7[15], stackBuffers[3].bundle(13));
		returnBuffer.bundle(13).multSub(valLocal6[14], stackBuffers[1].bundle(9));
		returnBuffer.bundle(13).multSub(valLocal7[14], stackBuffers[3].bundle(9));
		returnBuffer.bundle(13).multSub(valLocal6[13], stackBuffers[1].bundle(5));
		returnBuffer.bundle(13).multSub(valLocal7[13], stackBuffers[3].bundle(5));
		returnBuffer.bundle(13).multSub(valLocal6[12], stackBuffers[1].bundle(1));
		returnBuffer.bundle(13).multSub(valLocal7[12], stackBuffers[3].bundle(1));
		returnBuffer.bundle(13).multSub(valLocal6[11], stackBuffers[1].bundle(9));
		returnBuffer.bundle(13).multSub(valLocal7[11], stackBuffers[3].bundle(9));
		returnBuffer.bundle(13).multSub(valLocal6[10], stackBuffers[1].bundle(13));
		returnBuffer.bundle(13).multSub(valLocal7[10], stackBuffers[3].bundle(13));
		returnBuffer.bundle(13).multAdd(valLocal6[9], stackBuffers[1].bundle(1));
		returnBuffer.bundle(13).multAdd(valLocal7[9], stackBuffers[3].bundle(1));
		returnBuffer.bundle(13).multSub(valLocal6[8], stackBuffers[1].bundle(5));
		returnBuffer.bundle(13).multSub(valLocal7[8], stackBuffers[3].bundle(5));
		returnBuffer.bundle(13).multSub(valLocal6[7], stackBuffers[1].bundle(5));
		returnBuffer.bundle(13).multSub(valLocal7[7], stackBuffers[3].bundle(5));
		returnBuffer.bundle(13).multSub(valLocal6[6], stackBuffers[1].bundle(1));
		returnBuffer.bundle(13).multSub(valLocal7[6], stackBuffers[3].bundle(1));
		returnBuffer.bundle(13).multSub(valLocal6[5], stackBuffers[1].bundle(13));
		returnBuffer.bundle(13).multSub(valLocal7[5], stackBuffers[3].bundle(13));
		returnBuffer.bundle(13).multAdd(valLocal6[4], stackBuffers[1].bundle(9));
		returnBuffer.bundle(13).multAdd(valLocal7[4], stackBuffers[3].bundle(9));
		returnBuffer.bundle(13).multSub(valLocal6[3], stackBuffers[1].bundle(1));
		returnBuffer.bundle(13).multSub(valLocal7[3], stackBuffers[3].bundle(1));
		returnBuffer.bundle(13).multAdd(valLocal6[2], stackBuffers[1].bundle(5));
		returnBuffer.bundle(13).multAdd(valLocal7[2], stackBuffers[3].bundle(5));
		returnBuffer.bundle(13).multSub(valLocal6[1], stackBuffers[1].bundle(9));
		returnBuffer.bundle(13).multSub(valLocal7[1], stackBuffers[3].bundle(9));
		returnBuffer.bundle(13).multSub(valLocal6[0], stackBuffers[1].bundle(13));
		returnBuffer.bundle(13).multSub(valLocal7[0], stackBuffers[3].bundle(13));
		returnBuffer.bundle(14).multSub(valLocal6[15], stackBuffers[1].bundle(14));
		returnBuffer.bundle(14).multSub(valLocal7[15], stackBuffers[3].bundle(14));
		returnBuffer.bundle(14).multSub(valLocal6[14], stackBuffers[1].bundle(10));
		returnBuffer.bundle(14).multSub(valLocal7[14], stackBuffers[3].bundle(10));
		returnBuffer.bundle(14).multSub(valLocal6[13], stackBuffers[1].bundle(6));
		returnBuffer.bundle(14).multSub(valLocal7[13], stackBuffers[3].bundle(6));
		returnBuffer.bundle(14).multSub(valLocal6[12], stackBuffers[1].bundle(2));
		returnBuffer.bundle(14).multSub(valLocal7[12], stackBuffers[3].bundle(2));
		returnBuffer.bundle(14).multSub(valLocal6[11], stackBuffers[1].bundle(10));
		returnBuffer.bundle(14).multSub(valLocal7[11], stackBuffers[3].bundle(10));
		returnBuffer.bundle(14).multSub(valLocal6[10], stackBuffers[1].bundle(14));
		returnBuffer.bundle(14).multSub(valLocal7[10], stackBuffers[3].bundle(14));
		returnBuffer.bundle(14).multAdd(valLocal6[9], stackBuffers[1].bundle(2));
		returnBuffer.bundle(14).multAdd(valLocal7[9], stackBuffers[3].bundle(2));
		returnBuffer.bundle(14).multSub(valLocal6[8], stackBuffers[1].bundle(6));
		returnBuffer.bundle(14).multSub(valLocal7[8], stackBuffers[3].bundle(6));
		returnBuffer.bundle(14).multSub(valLocal6[7], stackBuffers[1].bundle(6));
		returnBuffer.bundle(14).multSub(valLocal7[7], stackBuffers[3].bundle(6));
		returnBuffer.bundle(14).multSub(valLocal6[6], stackBuffers[1].bundle(2));
		returnBuffer.bundle(14).multSub(valLocal7[6], stackBuffers[3].bundle(2));
		returnBuffer.bundle(14).multSub(valLocal6[5], stackBuffers[1].bundle(14));
		returnBuffer.bundle(14).multSub(valLocal7[5], stackBuffers[3].bundle(14));
		returnBuffer.bundle(14).multAdd(valLocal6[4], stackBuffers[1].bundle(10));
		returnBuffer.bundle(14).multAdd(valLocal7[4], stackBuffers[3].bundle(10));
		returnBuffer.bundle(14).multSub(valLocal6[3], stackBuffers[1].bundle(2));
		returnBuffer.bundle(14).multSub(valLocal7[3], stackBuffers[3].bundle(2));
		returnBuffer.bundle(14).multAdd(valLocal6[2], stackBuffers[1].bundle(6));
		returnBuffer.bundle(14).multAdd(valLocal7[2], stackBuffers[3].bundle(6));
		returnBuffer.bundle(14).multSub(valLocal6[1], stackBuffers[1].bundle(10));
		returnBuffer.bundle(14).multSub(valLocal7[1], stackBuffers[3].bundle(10));
		returnBuffer.bundle(14).multSub(valLocal6[0], stackBuffers[1].bundle(14));
		returnBuffer.bundle(14).multSub(valLocal7[0], stackBuffers[3].bundle(14));
		returnBuffer.bundle(3).multSub(valLocal6[15], stackBuffers[1].bundle(3));
		returnBuffer.bundle(3).multSub(valLocal7[15], stackBuffers[3].bundle(3));
		returnBuffer.bundle(3).multSub(valLocal6[14], stackBuffers[1].bundle(7));
		returnBuffer.bundle(3).multSub(valLocal7[14], stackBuffers[3].bundle(7));
		returnBuffer.bundle(3).multAdd(valLocal6[13], stackBuffers[1].bundle(11));
		returnBuffer.bundle(3).multAdd(valLocal7[13], stackBuffers[3].bundle(11));
		returnBuffer.bundle(3).multSub(valLocal6[12], stackBuffers[1].bundle(15));
		returnBuffer.bundle(3).multSub(valLocal7[12], stackBuffers[3].bundle(15));
		returnBuffer.bundle(3).multAdd(valLocal6[11], stackBuffers[1].bundle(7));
		returnBuffer.bundle(3).multAdd(valLocal7[11], stackBuffers[3].bundle(7));
		returnBuffer.bundle(3).multAdd(valLocal6[10], stackBuffers[1].bundle(3));
		returnBuffer.bundle(3).multAdd(valLocal7[10], stackBuffers[3].bundle(3));
		returnBuffer.bundle(3).multSub(valLocal6[9], stackBuffers[1].bundle(15));
		returnBuffer.bundle(3).multSub(valLocal7[9], stackBuffers[3].bundle(15));
		returnBuffer.bundle(3).multSub(valLocal6[8], stackBuffers[1].bundle(11));
		returnBuffer.bundle(3).multSub(valLocal7[8], stackBuffers[3].bundle(11));
		returnBuffer.bundle(3).multSub(valLocal6[7], stackBuffers[1].bundle(11));
		returnBuffer.bundle(3).multSub(valLocal7[7], stackBuffers[3].bundle(11));
		returnBuffer.bundle(3).multAdd(valLocal6[6], stackBuffers[1].bundle(15));
		returnBuffer.bundle(3).multAdd(valLocal7[6], stackBuffers[3].bundle(15));
		returnBuffer.bundle(3).multAdd(valLocal6[5], stackBuffers[1].bundle(3));
		returnBuffer.bundle(3).multAdd(valLocal7[5], stackBuffers[3].bundle(3));
		returnBuffer.bundle(3).multSub(valLocal6[4], stackBuffers[1].bundle(7));
		returnBuffer.bundle(3).multSub(valLocal7[4], stackBuffers[3].bundle(7));
		returnBuffer.bundle(3).multSub(valLocal6[3], stackBuffers[1].bundle(15));
		returnBuffer.bundle(3).multSub(valLocal7[3], stackBuffers[3].bundle(15));
		returnBuffer.bundle(3).multSub(valLocal6[2], stackBuffers[1].bundle(11));
		returnBuffer.bundle(3).multSub(valLocal7[2], stackBuffers[3].bundle(11));
		returnBuffer.bundle(3).multSub(valLocal6[1], stackBuffers[1].bundle(7));
		returnBuffer.bundle(3).multSub(valLocal7[1], stackBuffers[3].bundle(7));
		returnBuffer.bundle(3).multSub(valLocal6[0], stackBuffers[1].bundle(3));
		returnBuffer.bundle(3).multSub(valLocal7[0], stackBuffers[3].bundle(3));
		returnBuffer.bundle(0).multSub(valLocal6[15], stackBuffers[1].bundle(0));
		returnBuffer.bundle(0).multSub(valLocal7[15], stackBuffers[3].bundle(0));
		returnBuffer.bundle(0).multSub(valLocal6[14], stackBuffers[1].bundle(4));
		returnBuffer.bundle(0).multSub(valLocal7[14], stackBuffers[3].bundle(4));
		returnBuffer.bundle(0).multAdd(valLocal6[13], stackBuffers[1].bundle(8));
		returnBuffer.bundle(0).multAdd(valLocal7[13], stackBuffers[3].bundle(8));
		returnBuffer.bundle(0).multAdd(valLocal6[12], stackBuffers[1].bundle(12));
		returnBuffer.bundle(0).multAdd(valLocal7[12], stackBuffers[3].bundle(12));
		returnBuffer.bundle(0).multAdd(valLocal6[11], stackBuffers[1].bundle(4));
		returnBuffer.bundle(0).multAdd(valLocal7[11], stackBuffers[3].bundle(4));
		returnBuffer.bundle(0).multAdd(valLocal6[10], stackBuffers[1].bundle(0));
		returnBuffer.bundle(0).multAdd(valLocal7[10], stackBuffers[3].bundle(0));
		returnBuffer.bundle(0).multAdd(valLocal6[9], stackBuffers[1].bundle(12));
		returnBuffer.bundle(0).multAdd(valLocal7[9], stackBuffers[3].bundle(12));
		returnBuffer.bundle(0).multSub(valLocal6[8], stackBuffers[1].bundle(8));
		returnBuffer.bundle(0).multSub(valLocal7[8], stackBuffers[3].bundle(8));
		returnBuffer.bundle(0).multSub(valLocal6[7], stackBuffers[1].bundle(8));
		returnBuffer.bundle(0).multSub(valLocal7[7], stackBuffers[3].bundle(8));
		returnBuffer.bundle(0).multSub(valLocal6[6], stackBuffers[1].bundle(12));
		returnBuffer.bundle(0).multSub(valLocal7[6], stackBuffers[3].bundle(12));
		returnBuffer.bundle(0).multAdd(valLocal6[5], stackBuffers[1].bundle(0));
		returnBuffer.bundle(0).multAdd(valLocal7[5], stackBuffers[3].bundle(0));
		returnBuffer.bundle(0).multSub(valLocal6[4], stackBuffers[1].bundle(4));
		returnBuffer.bundle(0).multSub(valLocal7[4], stackBuffers[3].bundle(4));
		returnBuffer.bundle(0).multAdd(valLocal6[3], stackBuffers[1].bundle(12));
		returnBuffer.bundle(0).multAdd(valLocal7[3], stackBuffers[3].bundle(12));
		returnBuffer.bundle(0).multSub(valLocal6[2], stackBuffers[1].bundle(8));
		returnBuffer.bundle(0).multSub(valLocal7[2], stackBuffers[3].bundle(8));
		returnBuffer.bundle(0).multSub(valLocal6[1], stackBuffers[1].bundle(4));
		returnBuffer.bundle(0).multSub(valLocal7[1], stackBuffers[3].bundle(4));
		returnBuffer.bundle(0).multSub(valLocal6[0], stackBuffers[1].bundle(0));
		returnBuffer.bundle(0).multSub(valLocal7[0], stackBuffers[3].bundle(0));
		returnBuffer.bundle(1).multSub(valLocal6[15], stackBuffers[1].bundle(1));
		returnBuffer.bundle(1).multSub(valLocal7[15], stackBuffers[3].bundle(1));
		returnBuffer.bundle(1).multSub(valLocal6[14], stackBuffers[1].bundle(5));
		returnBuffer.bundle(1).multSub(valLocal7[14], stackBuffers[3].bundle(5));
		returnBuffer.bundle(1).multAdd(valLocal6[13], stackBuffers[1].bundle(9));
		returnBuffer.bundle(1).multAdd(valLocal7[13], stackBuffers[3].bundle(9));
		returnBuffer.bundle(1).multAdd(valLocal6[12], stackBuffers[1].bundle(13));
		returnBuffer.bundle(1).multAdd(valLocal7[12], stackBuffers[3].bundle(13));
		returnBuffer.bundle(1).multAdd(valLocal6[11], stackBuffers[1].bundle(5));
		returnBuffer.bundle(1).multAdd(valLocal7[11], stackBuffers[3].bundle(5));
		returnBuffer.bundle(1).multAdd(valLocal6[10], stackBuffers[1].bundle(1));
		returnBuffer.bundle(1).multAdd(valLocal7[10], stackBuffers[3].bundle(1));
		returnBuffer.bundle(1).multAdd(valLocal6[9], stackBuffers[1].bundle(13));
		returnBuffer.bundle(1).multAdd(valLocal7[9], stackBuffers[3].bundle(13));
		returnBuffer.bundle(1).multSub(valLocal6[8], stackBuffers[1].bundle(9));
		returnBuffer.bundle(1).multSub(valLocal7[8], stackBuffers[3].bundle(9));
		returnBuffer.bundle(1).multSub(valLocal6[7], stackBuffers[1].bundle(9));
		returnBuffer.bundle(1).multSub(valLocal7[7], stackBuffers[3].bundle(9));
		returnBuffer.bundle(1).multSub(valLocal6[6], stackBuffers[1].bundle(13));
		returnBuffer.bundle(1).multSub(valLocal7[6], stackBuffers[3].bundle(13));
		returnBuffer.bundle(1).multAdd(valLocal6[5], stackBuffers[1].bundle(1));
		returnBuffer.bundle(1).multAdd(valLocal7[5], stackBuffers[3].bundle(1));
		returnBuffer.bundle(1).multSub(valLocal6[4], stackBuffers[1].bundle(5));
		returnBuffer.bundle(1).multSub(valLocal7[4], stackBuffers[3].bundle(5));
		returnBuffer.bundle(1).multAdd(valLocal6[3], stackBuffers[1].bundle(13));
		returnBuffer.bundle(1).multAdd(valLocal7[3], stackBuffers[3].bundle(13));
		returnBuffer.bundle(1).multSub(valLocal6[2], stackBuffers[1].bundle(9));
		returnBuffer.bundle(1).multSub(valLocal7[2], stackBuffers[3].bundle(9));
		returnBuffer.bundle(1).multSub(valLocal6[1], stackBuffers[1].bundle(5));
		returnBuffer.bundle(1).multSub(valLocal7[1], stackBuffers[3].bundle(5));
		returnBuffer.bundle(1).multSub(valLocal6[0], stackBuffers[1].bundle(1));
		returnBuffer.bundle(1).multSub(valLocal7[0], stackBuffers[3].bundle(1));
		returnBuffer.bundle(2).multSub(valLocal6[15], stackBuffers[1].bundle(2));
		returnBuffer.bundle(2).multSub(valLocal7[15], stackBuffers[3].bundle(2));
		returnBuffer.bundle(2).multSub(valLocal6[14], stackBuffers[1].bundle(6));
		returnBuffer.bundle(2).multSub(valLocal7[14], stackBuffers[3].bundle(6));
		returnBuffer.bundle(2).multAdd(valLocal6[13], stackBuffers[1].bundle(10));
		returnBuffer.bundle(2).multAdd(valLocal7[13], stackBuffers[3].bundle(10));
		returnBuffer.bundle(2).multAdd(valLocal6[12], stackBuffers[1].bundle(14));
		returnBuffer.bundle(2).multAdd(valLocal7[12], stackBuffers[3].bundle(14));
		returnBuffer.bundle(2).multAdd(valLocal6[11], stackBuffers[1].bundle(6));
		returnBuffer.bundle(2).multAdd(valLocal7[11], stackBuffers[3].bundle(6));
		returnBuffer.bundle(2).multAdd(valLocal6[10], stackBuffers[1].bundle(2));
		returnBuffer.bundle(2).multAdd(valLocal7[10], stackBuffers[3].bundle(2));
		returnBuffer.bundle(2).multAdd(valLocal6[9], stackBuffers[1].bundle(14));
		returnBuffer.bundle(2).multAdd(valLocal7[9], stackBuffers[3].bundle(14));
		returnBuffer.bundle(2).multSub(valLocal6[8], stackBuffers[1].bundle(10));
		returnBuffer.bundle(2).multSub(valLocal7[8], stackBuffers[3].bundle(10));
		returnBuffer.bundle(2).multSub(valLocal6[7], stackBuffers[1].bundle(10));
		returnBuffer.bundle(2).multSub(valLocal7[7], stackBuffers[3].bundle(10));
		returnBuffer.bundle(2).multSub(valLocal6[6], stackBuffers[1].bundle(14));
		returnBuffer.bundle(2).multSub(valLocal7[6], stackBuffers[3].bundle(14));
		returnBuffer.bundle(2).multAdd(valLocal6[5], stackBuffers[1].bundle(2));
		returnBuffer.bundle(2).multAdd(valLocal7[5], stackBuffers[3].bundle(2));
		returnBuffer.bundle(2).multSub(valLocal6[4], stackBuffers[1].bundle(6));
		returnBuffer.bundle(2).multSub(valLocal7[4], stackBuffers[3].bundle(6));
		returnBuffer.bundle(2).multAdd(valLocal6[3], stackBuffers[1].bundle(14));
		returnBuffer.bundle(2).multAdd(valLocal7[3], stackBuffers[3].bundle(14));
		returnBuffer.bundle(2).multSub(valLocal6[2], stackBuffers[1].bundle(10));
		returnBuffer.bundle(2).multSub(valLocal7[2], stackBuffers[3].bundle(10));
		returnBuffer.bundle(2).multSub(valLocal6[1], stackBuffers[1].bundle(6));
		returnBuffer.bundle(2).multSub(valLocal7[1], stackBuffers[3].bundle(6));
		returnBuffer.bundle(2).multSub(valLocal6[0], stackBuffers[1].bundle(2));
		returnBuffer.bundle(2).multSub(valLocal7[0], stackBuffers[3].bundle(2));
		returnBuffer.bundle(7).multSub(valLocal6[15], stackBuffers[1].bundle(7));
		returnBuffer.bundle(7).multSub(valLocal7[15], stackBuffers[3].bundle(7));
		returnBuffer.bundle(7).multAdd(valLocal6[14], stackBuffers[1].bundle(3));
		returnBuffer.bundle(7).multAdd(valLocal7[14], stackBuffers[3].bundle(3));
		returnBuffer.bundle(7).multSub(valLocal6[13], stackBuffers[1].bundle(15));
		returnBuffer.bundle(7).multSub(valLocal7[13], stackBuffers[3].bundle(15));
		returnBuffer.bundle(7).multSub(valLocal6[12], stackBuffers[1].bundle(11));
		returnBuffer.bundle(7).multSub(valLocal7[12], stackBuffers[3].bundle(11));
		returnBuffer.bundle(7).multSub(valLocal6[11], stackBuffers[1].bundle(3));
		returnBuffer.bundle(7).multSub(valLocal7[11], stackBuffers[3].bundle(3));
		returnBuffer.bundle(7).multAdd(valLocal6[10], stackBuffers[1].bundle(7));
		returnBuffer.bundle(7).multAdd(valLocal7[10], stackBuffers[3].bundle(7));
		returnBuffer.bundle(7).multSub(valLocal6[9], stackBuffers[1].bundle(11));
		returnBuffer.bundle(7).multSub(valLocal7[9], stackBuffers[3].bundle(11));
		returnBuffer.bundle(7).multAdd(valLocal6[8], stackBuffers[1].bundle(15));
		returnBuffer.bundle(7).multAdd(valLocal7[8], stackBuffers[3].bundle(15));
		returnBuffer.bundle(7).multSub(valLocal6[7], stackBuffers[1].bundle(15));
		returnBuffer.bundle(7).multSub(valLocal7[7], stackBuffers[3].bundle(15));
		returnBuffer.bundle(7).multSub(valLocal6[6], stackBuffers[1].bundle(11));
		returnBuffer.bundle(7).multSub(valLocal7[6], stackBuffers[3].bundle(11));
		returnBuffer.bundle(7).multSub(valLocal6[5], stackBuffers[1].bundle(7));
		returnBuffer.bundle(7).multSub(valLocal7[5], stackBuffers[3].bundle(7));
		returnBuffer.bundle(7).multSub(valLocal6[4], stackBuffers[1].bundle(3));
		returnBuffer.bundle(7).multSub(valLocal7[4], stackBuffers[3].bundle(3));
		returnBuffer.bundle(7).multAdd(valLocal6[3], stackBuffers[1].bundle(11));
		returnBuffer.bundle(7).multAdd(valLocal7[3], stackBuffers[3].bundle(11));
		returnBuffer.bundle(7).multSub(valLocal6[2], stackBuffers[1].bundle(15));
		returnBuffer.bundle(7).multSub(valLocal7[2], stackBuffers[3].bundle(15));
		returnBuffer.bundle(7).multSub(valLocal6[1], stackBuffers[1].bundle(3));
		returnBuffer.bundle(7).multSub(valLocal7[1], stackBuffers[3].bundle(3));
		returnBuffer.bundle(7).multAdd(valLocal6[0], stackBuffers[1].bundle(7));
		returnBuffer.bundle(7).multAdd(valLocal7[0], stackBuffers[3].bundle(7));
		returnBuffer.bundle(4).multSub(valLocal6[15], stackBuffers[1].bundle(4));
		returnBuffer.bundle(4).multSub(valLocal7[15], stackBuffers[3].bundle(4));
		returnBuffer.bundle(4).multAdd(valLocal6[14], stackBuffers[1].bundle(0));
		returnBuffer.bundle(4).multAdd(valLocal7[14], stackBuffers[3].bundle(0));
		returnBuffer.bundle(4).multAdd(valLocal6[13], stackBuffers[1].bundle(12));
		returnBuffer.bundle(4).multAdd(valLocal7[13], stackBuffers[3].bundle(12));
		returnBuffer.bundle(4).multSub(valLocal6[12], stackBuffers[1].bundle(8));
		returnBuffer.bundle(4).multSub(valLocal7[12], stackBuffers[3].bundle(8));
		returnBuffer.bundle(4).multSub(valLocal6[11], stackBuffers[1].bundle(0));
		returnBuffer.bundle(4).multSub(valLocal7[11], stackBuffers[3].bundle(0));
		returnBuffer.bundle(4).multAdd(valLocal6[10], stackBuffers[1].bundle(4));
		returnBuffer.bundle(4).multAdd(valLocal7[10], stackBuffers[3].bundle(4));
		returnBuffer.bundle(4).multSub(valLocal6[9], stackBuffers[1].bundle(8));
		returnBuffer.bundle(4).multSub(valLocal7[9], stackBuffers[3].bundle(8));
		returnBuffer.bundle(4).multSub(valLocal6[8], stackBuffers[1].bundle(12));
		returnBuffer.bundle(4).multSub(valLocal7[8], stackBuffers[3].bundle(12));
		returnBuffer.bundle(4).multAdd(valLocal6[7], stackBuffers[1].bundle(12));
		returnBuffer.bundle(4).multAdd(valLocal7[7], stackBuffers[3].bundle(12));
		returnBuffer.bundle(4).multSub(valLocal6[6], stackBuffers[1].bundle(8));
		returnBuffer.bundle(4).multSub(valLocal7[6], stackBuffers[3].bundle(8));
		returnBuffer.bundle(4).multSub(valLocal6[5], stackBuffers[1].bundle(4));
		returnBuffer.bundle(4).multSub(valLocal7[5], stackBuffers[3].bundle(4));
		returnBuffer.bundle(4).multSub(valLocal6[4], stackBuffers[1].bundle(0));
		returnBuffer.bundle(4).multSub(valLocal7[4], stackBuffers[3].bundle(0));
		returnBuffer.bundle(4).multAdd(valLocal6[3], stackBuffers[1].bundle(8));
		returnBuffer.bundle(4).multAdd(valLocal7[3], stackBuffers[3].bundle(8));
		returnBuffer.bundle(4).multAdd(valLocal6[2], stackBuffers[1].bundle(12));
		returnBuffer.bundle(4).multAdd(valLocal7[2], stackBuffers[3].bundle(12));
		returnBuffer.bundle(4).multSub(valLocal6[1], stackBuffers[1].bundle(0));
		returnBuffer.bundle(4).multSub(valLocal7[1], stackBuffers[3].bundle(0));
		returnBuffer.bundle(4).multAdd(valLocal6[0], stackBuffers[1].bundle(4));
		returnBuffer.bundle(4).multAdd(valLocal7[0], stackBuffers[3].bundle(4));
		returnBuffer.bundle(5).multSub(valLocal6[15], stackBuffers[1].bundle(5));
		returnBuffer.bundle(5).multSub(valLocal7[15], stackBuffers[3].bundle(5));
		returnBuffer.bundle(5).multAdd(valLocal6[14], stackBuffers[1].bundle(1));
		returnBuffer.bundle(5).multAdd(valLocal7[14], stackBuffers[3].bundle(1));
		returnBuffer.bundle(5).multAdd(valLocal6[13], stackBuffers[1].bundle(13));
		returnBuffer.bundle(5).multAdd(valLocal7[13], stackBuffers[3].bundle(13));
		returnBuffer.bundle(5).multSub(valLocal6[12], stackBuffers[1].bundle(9));
		returnBuffer.bundle(5).multSub(valLocal7[12], stackBuffers[3].bundle(9));
		returnBuffer.bundle(5).multSub(valLocal6[11], stackBuffers[1].bundle(1));
		returnBuffer.bundle(5).multSub(valLocal7[11], stackBuffers[3].bundle(1));
		returnBuffer.bundle(5).multAdd(valLocal6[10], stackBuffers[1].bundle(5));
		returnBuffer.bundle(5).multAdd(valLocal7[10], stackBuffers[3].bundle(5));
		returnBuffer.bundle(5).multSub(valLocal6[9], stackBuffers[1].bundle(9));
		returnBuffer.bundle(5).multSub(valLocal7[9], stackBuffers[3].bundle(9));
		returnBuffer.bundle(5).multSub(valLocal6[8], stackBuffers[1].bundle(13));
		returnBuffer.bundle(5).multSub(valLocal7[8], stackBuffers[3].bundle(13));
		returnBuffer.bundle(5).multAdd(valLocal6[7], stackBuffers[1].bundle(13));
		returnBuffer.bundle(5).multAdd(valLocal7[7], stackBuffers[3].bundle(13));
		returnBuffer.bundle(5).multSub(valLocal6[6], stackBuffers[1].bundle(9));
		returnBuffer.bundle(5).multSub(valLocal7[6], stackBuffers[3].bundle(9));
		returnBuffer.bundle(5).multSub(valLocal6[5], stackBuffers[1].bundle(5));
		returnBuffer.bundle(5).multSub(valLocal7[5], stackBuffers[3].bundle(5));
		returnBuffer.bundle(5).multSub(valLocal6[4], stackBuffers[1].bundle(1));
		returnBuffer.bundle(5).multSub(valLocal7[4], stackBuffers[3].bundle(1));
		returnBuffer.bundle(5).multAdd(valLocal6[3], stackBuffers[1].bundle(9));
		returnBuffer.bundle(5).multAdd(valLocal7[3], stackBuffers[3].bundle(9));
		returnBuffer.bundle(5).multAdd(valLocal6[2], stackBuffers[1].bundle(13));
		returnBuffer.bundle(5).multAdd(valLocal7[2], stackBuffers[3].bundle(13));
		returnBuffer.bundle(5).multSub(valLocal6[1], stackBuffers[1].bundle(1));
		returnBuffer.bundle(5).multSub(valLocal7[1], stackBuffers[3].bundle(1));
		returnBuffer.bundle(5).multAdd(valLocal6[0], stackBuffers[1].bundle(5));
		returnBuffer.bundle(5).multAdd(valLocal7[0], stackBuffers[3].bundle(5));
		returnBuffer.bundle(6).multSub(valLocal6[15], stackBuffers[1].bundle(6));
		returnBuffer.bundle(6).multSub(valLocal7[15], stackBuffers[3].bundle(6));
		returnBuffer.bundle(6).multAdd(valLocal6[14], stackBuffers[1].bundle(2));
		returnBuffer.bundle(6).multAdd(valLocal7[14], stackBuffers[3].bundle(2));
		returnBuffer.bundle(6).multAdd(valLocal6[13], stackBuffers[1].bundle(14));
		returnBuffer.bundle(6).multAdd(valLocal7[13], stackBuffers[3].bundle(14));
		returnBuffer.bundle(6).multSub(valLocal6[12], stackBuffers[1].bundle(10));
		returnBuffer.bundle(6).multSub(valLocal7[12], stackBuffers[3].bundle(10));
		returnBuffer.bundle(6).multSub(valLocal6[11], stackBuffers[1].bundle(2));
		returnBuffer.bundle(6).multSub(valLocal7[11], stackBuffers[3].bundle(2));
		returnBuffer.bundle(6).multAdd(valLocal6[10], stackBuffers[1].bundle(6));
		returnBuffer.bundle(6).multAdd(valLocal7[10], stackBuffers[3].bundle(6));
		returnBuffer.bundle(6).multSub(valLocal6[9], stackBuffers[1].bundle(10));
		returnBuffer.bundle(6).multSub(valLocal7[9], stackBuffers[3].bundle(10));
		returnBuffer.bundle(6).multSub(valLocal6[8], stackBuffers[1].bundle(14));
		returnBuffer.bundle(6).multSub(valLocal7[8], stackBuffers[3].bundle(14));
		returnBuffer.bundle(6).multAdd(valLocal6[7], stackBuffers[1].bundle(14));
		returnBuffer.bundle(6).multAdd(valLocal7[7], stackBuffers[3].bundle(14));
		returnBuffer.bundle(6).multSub(valLocal6[6], stackBuffers[1].bundle(10));
		returnBuffer.bundle(6).multSub(valLocal7[6], stackBuffers[3].bundle(10));
		returnBuffer.bundle(6).multSub(valLocal6[5], stackBuffers[1].bundle(6));
		returnBuffer.bundle(6).multSub(valLocal7[5], stackBuffers[3].bundle(6));
		returnBuffer.bundle(6).multSub(valLocal6[4], stackBuffers[1].bundle(2));
		returnBuffer.bundle(6).multSub(valLocal7[4], stackBuffers[3].bundle(2));
		returnBuffer.bundle(6).multAdd(valLocal6[3], stackBuffers[1].bundle(10));
		returnBuffer.bundle(6).multAdd(valLocal7[3], stackBuffers[3].bundle(10));
		returnBuffer.bundle(6).multAdd(valLocal6[2], stackBuffers[1].bundle(14));
		returnBuffer.bundle(6).multAdd(valLocal7[2], stackBuffers[3].bundle(14));
		returnBuffer.bundle(6).multSub(valLocal6[1], stackBuffers[1].bundle(2));
		returnBuffer.bundle(6).multSub(valLocal7[1], stackBuffers[3].bundle(2));
		returnBuffer.bundle(6).multAdd(valLocal6[0], stackBuffers[1].bundle(6));
		returnBuffer.bundle(6).multAdd(valLocal7[0], stackBuffers[3].bundle(6));
		returnBuffer.bundle(11).multSub(valLocal6[15], stackBuffers[1].bundle(11));
		returnBuffer.bundle(11).multSub(valLocal7[15], stackBuffers[3].bundle(11));
		returnBuffer.bundle(11).multSub(valLocal6[14], stackBuffers[1].bundle(15));
		returnBuffer.bundle(11).multSub(valLocal7[14], stackBuffers[3].bundle(15));
		returnBuffer.bundle(11).multSub(valLocal6[13], stackBuffers[1].bundle(3));
		returnBuffer.bundle(11).multSub(valLocal7[13], stackBuffers[3].bundle(3));
		returnBuffer.bundle(11).multAdd(valLocal6[12], stackBuffers[1].bundle(7));
		returnBuffer.bundle(11).multAdd(valLocal7[12], stackBuffers[3].bundle(7));
		returnBuffer.bundle(11).multSub(valLocal6[11], stackBuffers[1].bundle(15));
		returnBuffer.bundle(11).multSub(valLocal7[11], stackBuffers[3].bundle(15));
		returnBuffer.bundle(11).multSub(valLocal6[10], stackBuffers[1].bundle(11));
		returnBuffer.bundle(11).multSub(valLocal7[10], stackBuffers[3].bundle(11));
		returnBuffer.bundle(11).multSub(valLocal6[9], stackBuffers[1].bundle(7));
		returnBuffer.bundle(11).multSub(valLocal7[9], stackBuffers[3].bundle(7));
		returnBuffer.bundle(11).multSub(valLocal6[8], stackBuffers[1].bundle(3));
		returnBuffer.bundle(11).multSub(valLocal7[8], stackBuffers[3].bundle(3));
		returnBuffer.bundle(11).multAdd(valLocal6[7], stackBuffers[1].bundle(3));
		returnBuffer.bundle(11).multAdd(valLocal7[7], stackBuffers[3].bundle(3));
		returnBuffer.bundle(11).multSub(valLocal6[6], stackBuffers[1].bundle(7));
		returnBuffer.bundle(11).multSub(valLocal7[6], stackBuffers[3].bundle(7));
		returnBuffer.bundle(11).multAdd(valLocal6[5], stackBuffers[1].bundle(11));
		returnBuffer.bundle(11).multAdd(valLocal7[5], stackBuffers[3].bundle(11));
		returnBuffer.bundle(11).multSub(valLocal6[4], stackBuffers[1].bundle(15));
		returnBuffer.bundle(11).multSub(valLocal7[4], stackBuffers[3].bundle(15));
		returnBuffer.bundle(11).multSub(valLocal6[3], stackBuffers[1].bundle(7));
		returnBuffer.bundle(11).multSub(valLocal7[3], stackBuffers[3].bundle(7));
		returnBuffer.bundle(11).multSub(valLocal6[2], stackBuffers[1].bundle(3));
		returnBuffer.bundle(11).multSub(valLocal7[2], stackBuffers[3].bundle(3));
		returnBuffer.bundle(11).multAdd(valLocal6[1], stackBuffers[1].bundle(15));
		returnBuffer.bundle(11).multAdd(valLocal7[1], stackBuffers[3].bundle(15));
		returnBuffer.bundle(11).multAdd(valLocal6[0], stackBuffers[1].bundle(11));
		returnBuffer.bundle(11).multAdd(valLocal7[0], stackBuffers[3].bundle(11));
		returnBuffer.bundle(8).multSub(valLocal6[15], stackBuffers[1].bundle(8));
		returnBuffer.bundle(8).multSub(valLocal7[15], stackBuffers[3].bundle(8));
		returnBuffer.bundle(8).multAdd(valLocal6[14], stackBuffers[1].bundle(12));
		returnBuffer.bundle(8).multAdd(valLocal7[14], stackBuffers[3].bundle(12));
		returnBuffer.bundle(8).multSub(valLocal6[13], stackBuffers[1].bundle(0));
		returnBuffer.bundle(8).multSub(valLocal7[13], stackBuffers[3].bundle(0));
		returnBuffer.bundle(8).multAdd(valLocal6[12], stackBuffers[1].bundle(4));
		returnBuffer.bundle(8).multAdd(valLocal7[12], stackBuffers[3].bundle(4));
		returnBuffer.bundle(8).multAdd(valLocal6[11], stackBuffers[1].bundle(12));
		returnBuffer.bundle(8).multAdd(valLocal7[11], stackBuffers[3].bundle(12));
		returnBuffer.bundle(8).multSub(valLocal6[10], stackBuffers[1].bundle(8));
		returnBuffer.bundle(8).multSub(valLocal7[10], stackBuffers[3].bundle(8));
		returnBuffer.bundle(8).multSub(valLocal6[9], stackBuffers[1].bundle(4));
		returnBuffer.bundle(8).multSub(valLocal7[9], stackBuffers[3].bundle(4));
		returnBuffer.bundle(8).multSub(valLocal6[8], stackBuffers[1].bundle(0));
		returnBuffer.bundle(8).multSub(valLocal7[8], stackBuffers[3].bundle(0));
		returnBuffer.bundle(8).multAdd(valLocal6[7], stackBuffers[1].bundle(0));
		returnBuffer.bundle(8).multAdd(valLocal7[7], stackBuffers[3].bundle(0));
		returnBuffer.bundle(8).multSub(valLocal6[6], stackBuffers[1].bundle(4));
		returnBuffer.bundle(8).multSub(valLocal7[6], stackBuffers[3].bundle(4));
		returnBuffer.bundle(8).multAdd(valLocal6[5], stackBuffers[1].bundle(8));
		returnBuffer.bundle(8).multAdd(valLocal7[5], stackBuffers[3].bundle(8));
		returnBuffer.bundle(8).multAdd(valLocal6[4], stackBuffers[1].bundle(12));
		returnBuffer.bundle(8).multAdd(valLocal7[4], stackBuffers[3].bundle(12));
		returnBuffer.bundle(8).multSub(valLocal6[3], stackBuffers[1].bundle(4));
		returnBuffer.bundle(8).multSub(valLocal7[3], stackBuffers[3].bundle(4));
		returnBuffer.bundle(8).multSub(valLocal6[2], stackBuffers[1].bundle(0));
		returnBuffer.bundle(8).multSub(valLocal7[2], stackBuffers[3].bundle(0));
		returnBuffer.bundle(8).multSub(valLocal6[1], stackBuffers[1].bundle(12));
		returnBuffer.bundle(8).multSub(valLocal7[1], stackBuffers[3].bundle(12));
		returnBuffer.bundle(8).multAdd(valLocal6[0], stackBuffers[1].bundle(8));
		returnBuffer.bundle(8).multAdd(valLocal7[0], stackBuffers[3].bundle(8));
		returnBuffer.bundle(9).multSub(valLocal6[15], stackBuffers[1].bundle(9));
		returnBuffer.bundle(9).multSub(valLocal7[15], stackBuffers[3].bundle(9));
		returnBuffer.bundle(9).multAdd(valLocal6[14], stackBuffers[1].bundle(13));
		returnBuffer.bundle(9).multAdd(valLocal7[14], stackBuffers[3].bundle(13));
		returnBuffer.bundle(9).multSub(valLocal6[13], stackBuffers[1].bundle(1));
		returnBuffer.bundle(9).multSub(valLocal7[13], stackBuffers[3].bundle(1));
		returnBuffer.bundle(9).multAdd(valLocal6[12], stackBuffers[1].bundle(5));
		returnBuffer.bundle(9).multAdd(valLocal7[12], stackBuffers[3].bundle(5));
		returnBuffer.bundle(9).multAdd(valLocal6[11], stackBuffers[1].bundle(13));
		returnBuffer.bundle(9).multAdd(valLocal7[11], stackBuffers[3].bundle(13));
		returnBuffer.bundle(9).multSub(valLocal6[10], stackBuffers[1].bundle(9));
		returnBuffer.bundle(9).multSub(valLocal7[10], stackBuffers[3].bundle(9));
		returnBuffer.bundle(9).multSub(valLocal6[9], stackBuffers[1].bundle(5));
		returnBuffer.bundle(9).multSub(valLocal7[9], stackBuffers[3].bundle(5));
		returnBuffer.bundle(9).multSub(valLocal6[8], stackBuffers[1].bundle(1));
		returnBuffer.bundle(9).multSub(valLocal7[8], stackBuffers[3].bundle(1));
		returnBuffer.bundle(9).multAdd(valLocal6[7], stackBuffers[1].bundle(1));
		returnBuffer.bundle(9).multAdd(valLocal7[7], stackBuffers[3].bundle(1));
		returnBuffer.bundle(9).multSub(valLocal6[6], stackBuffers[1].bundle(5));
		returnBuffer.bundle(9).multSub(valLocal7[6], stackBuffers[3].bundle(5));
		returnBuffer.bundle(9).multAdd(valLocal6[5], stackBuffers[1].bundle(9));
		returnBuffer.bundle(9).multAdd(valLocal7[5], stackBuffers[3].bundle(9));
		returnBuffer.bundle(9).multAdd(valLocal6[4], stackBuffers[1].bundle(13));
		returnBuffer.bundle(9).multAdd(valLocal7[4], stackBuffers[3].bundle(13));
		returnBuffer.bundle(9).multSub(valLocal6[3], stackBuffers[1].bundle(5));
		returnBuffer.bundle(9).multSub(valLocal7[3], stackBuffers[3].bundle(5));
		returnBuffer.bundle(9).multSub(valLocal6[2], stackBuffers[1].bundle(1));
		returnBuffer.bundle(9).multSub(valLocal7[2], stackBuffers[3].bundle(1));
		returnBuffer.bundle(9).multSub(valLocal6[1], stackBuffers[1].bundle(13));
		returnBuffer.bundle(9).multSub(valLocal7[1], stackBuffers[3].bundle(13));
		returnBuffer.bundle(9).multAdd(valLocal6[0], stackBuffers[1].bundle(9));
		returnBuffer.bundle(9).multAdd(valLocal7[0], stackBuffers[3].bundle(9));
		returnBuffer.bundle(10).multSub(valLocal6[15], stackBuffers[1].bundle(10));
		returnBuffer.bundle(10).multSub(valLocal7[15], stackBuffers[3].bundle(10));
		returnBuffer.bundle(10).multAdd(valLocal6[14], stackBuffers[1].bundle(14));
		returnBuffer.bundle(10).multAdd(valLocal7[14], stackBuffers[3].bundle(14));
		returnBuffer.bundle(10).multSub(valLocal6[13], stackBuffers[1].bundle(2));
		returnBuffer.bundle(10).multSub(valLocal7[13], stackBuffers[3].bundle(2));
		returnBuffer.bundle(10).multAdd(valLocal6[12], stackBuffers[1].bundle(6));
		returnBuffer.bundle(10).multAdd(valLocal7[12], stackBuffers[3].bundle(6));
		returnBuffer.bundle(10).multAdd(valLocal6[11], stackBuffers[1].bundle(14));
		returnBuffer.bundle(10).multAdd(valLocal7[11], stackBuffers[3].bundle(14));
		returnBuffer.bundle(10).multSub(valLocal6[10], stackBuffers[1].bundle(10));
		returnBuffer.bundle(10).multSub(valLocal7[10], stackBuffers[3].bundle(10));
		returnBuffer.bundle(10).multSub(valLocal6[9], stackBuffers[1].bundle(6));
		returnBuffer.bundle(10).multSub(valLocal7[9], stackBuffers[3].bundle(6));
		returnBuffer.bundle(10).multSub(valLocal6[8], stackBuffers[1].bundle(2));
		returnBuffer.bundle(10).multSub(valLocal7[8], stackBuffers[3].bundle(2));
		returnBuffer.bundle(10).multAdd(valLocal6[7], stackBuffers[1].bundle(2));
		returnBuffer.bundle(10).multAdd(valLocal7[7], stackBuffers[3].bundle(2));
		returnBuffer.bundle(10).multSub(valLocal6[6], stackBuffers[1].bundle(6));
		returnBuffer.bundle(10).multSub(valLocal7[6], stackBuffers[3].bundle(6));
		returnBuffer.bundle(10).multAdd(valLocal6[5], stackBuffers[1].bundle(10));
		returnBuffer.bundle(10).multAdd(valLocal7[5], stackBuffers[3].bundle(10));
		returnBuffer.bundle(10).multAdd(valLocal6[4], stackBuffers[1].bundle(14));
		returnBuffer.bundle(10).multAdd(valLocal7[4], stackBuffers[3].bundle(14));
		returnBuffer.bundle(10).multSub(valLocal6[3], stackBuffers[1].bundle(6));
		returnBuffer.bundle(10).multSub(valLocal7[3], stackBuffers[3].bundle(6));
		returnBuffer.bundle(10).multSub(valLocal6[2], stackBuffers[1].bundle(2));
		returnBuffer.bundle(10).multSub(valLocal7[2], stackBuffers[3].bundle(2));
		returnBuffer.bundle(10).multSub(valLocal6[1], stackBuffers[1].bundle(14));
		returnBuffer.bundle(10).multSub(valLocal7[1], stackBuffers[3].bundle(14));
		returnBuffer.bundle(10).multAdd(valLocal6[0], stackBuffers[1].bundle(10));
		returnBuffer.bundle(10).multAdd(valLocal7[0], stackBuffers[3].bundle(10));
		#pragma endregion
	};

	auto integralKernelU = [&](const float wp, ValueSuperbundle<float, 16> &returnBuffer) -> void
	{
		//u-Channel, to be combined with P(wp, u + wp) + P(u + wp, wp)
		//ph-ladder A and B, respectively (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab0 = v4->generateAccessBuffer(w1 + wp, -w2p + wp, u, TRIVertexTwoParticle::FrequencyChannel::U);
		const TRIVertexTwoParticleAccessBuffer<4> ab1 = v4->generateAccessBuffer(w1p + wp, w2 - wp, u, TRIVertexTwoParticle::FrequencyChannel::U);
		//ph-ladder A and B, respectively (negative sign)
		const TRIVertexTwoParticleAccessBuffer<4> ab2 = v4->generateAccessBuffer(w2p - wp, -w1 - wp, u, TRIVertexTwoParticle::FrequencyChannel::U);
		const TRIVertexTwoParticleAccessBuffer<4> ab3 = v4->generateAccessBuffer(w2 - wp, w1p + wp, u, TRIVertexTwoParticle::FrequencyChannel::U);

		v4->getValueSuperbundle(ab0, stackBuffers[0]);
		v4->getValueSuperbundle(ab1, stackBuffers[1]);
		v4->getValueSuperbundle(ab2, stackBuffers[2]);
		v4->getValueSuperbundle(ab3, stackBuffers[3]);

		//calculate _flow
		returnBuffer.reset();

		#pragma region phLadder
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(15));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(15));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(14));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(14));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(13));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(13));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(12));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(12));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(11));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(11));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(10));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(10));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(9));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(9));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(8));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(8));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(7));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(7));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(6));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(6));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(5));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(5));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(4));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(4));
		returnBuffer.bundle(15).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(3));
		returnBuffer.bundle(15).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(3));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(2));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(2));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(1));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(1));
		returnBuffer.bundle(15).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(0));
		returnBuffer.bundle(15).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(0));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(12));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(12));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(13));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(13));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(14));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(14));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(15));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(15));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(8));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(8));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(9));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(9));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(10));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(10));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(11));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(11));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(4));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(4));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(5));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(5));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(6));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(6));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(7));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(7));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(0));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(0));
		returnBuffer.bundle(12).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(1));
		returnBuffer.bundle(12).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(1));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(2));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(2));
		returnBuffer.bundle(12).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(3));
		returnBuffer.bundle(12).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(3));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(13));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(13));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(12));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(12));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(15));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(15));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(14));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(14));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(9));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(9));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(8));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(8));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(11));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(11));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(10));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(10));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(5));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(5));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(4));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(4));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(7));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(7));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(6));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(6));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(1));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(1));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(0));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(0));
		returnBuffer.bundle(13).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(3));
		returnBuffer.bundle(13).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(3));
		returnBuffer.bundle(13).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(2));
		returnBuffer.bundle(13).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(2));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(14));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(14));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(15));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(15));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(12));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(12));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(13));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(13));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(10));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(10));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(11));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(11));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(8));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(8));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(9));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(9));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(6));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(6));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(7));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(7));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(4));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(4));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(5));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(5));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(2));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(2));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(3));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(3));
		returnBuffer.bundle(14).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(0));
		returnBuffer.bundle(14).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(0));
		returnBuffer.bundle(14).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(1));
		returnBuffer.bundle(14).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(1));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(3));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(3));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(2));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(2));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(1));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(1));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(0));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(0));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(7));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(7));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(6));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(6));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(5));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(5));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(4));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(4));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(11));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(11));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(10));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(10));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(9));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(9));
		returnBuffer.bundle(3).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(8));
		returnBuffer.bundle(3).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(8));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(15));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(15));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(14));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(14));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(13));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(13));
		returnBuffer.bundle(3).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(12));
		returnBuffer.bundle(3).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(12));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(0));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(0));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(1));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(1));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(2));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(2));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(3));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(3));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(4));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(4));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(5));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(5));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(6));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(6));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(7));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(7));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(8));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(8));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(9));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(9));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(10));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(10));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(11));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(11));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(12));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(12));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(13));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(13));
		returnBuffer.bundle(0).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(14));
		returnBuffer.bundle(0).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(14));
		returnBuffer.bundle(0).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(15));
		returnBuffer.bundle(0).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(15));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(1));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(1));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(0));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(0));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(3));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(3));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(2));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(2));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(5));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(5));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(4));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(4));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(7));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(7));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(6));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(6));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(9));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(9));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(8));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(8));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(11));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(11));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(10));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(10));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(13));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(13));
		returnBuffer.bundle(1).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(12));
		returnBuffer.bundle(1).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(12));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(15));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(15));
		returnBuffer.bundle(1).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(14));
		returnBuffer.bundle(1).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(14));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(2));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(2));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(3));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(3));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(0));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(0));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(1));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(1));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(6));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(6));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(7));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(7));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(4));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(4));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(5));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(5));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(10));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(10));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(11));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(11));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(8));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(8));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(9));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(9));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(14));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(14));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(15));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(15));
		returnBuffer.bundle(2).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(12));
		returnBuffer.bundle(2).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(12));
		returnBuffer.bundle(2).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(13));
		returnBuffer.bundle(2).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(13));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(7));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(7));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(6));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(6));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(5));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(5));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(4));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(4));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(3));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(3));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(2));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(2));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(1));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(1));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(0));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(0));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(7), stackBuffers[1].bundle(15));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(7), stackBuffers[3].bundle(15));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(14));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(14));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(13));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(13));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(12));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(12));
		returnBuffer.bundle(7).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(11));
		returnBuffer.bundle(7).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(11));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(10));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(10));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(9));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(9));
		returnBuffer.bundle(7).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(8));
		returnBuffer.bundle(7).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(8));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(4));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(4));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(5));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(5));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(6));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(6));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(7));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(7));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(0));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(0));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(1));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(1));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(2));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(2));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(3));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(3));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(12));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(12));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(13));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(13));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(14));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(14));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(15));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(15));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(8));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(8));
		returnBuffer.bundle(4).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(9));
		returnBuffer.bundle(4).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(9));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(10));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(10));
		returnBuffer.bundle(4).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(11));
		returnBuffer.bundle(4).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(11));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(5));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(5));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(4));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(4));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(7));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(7));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(6));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(6));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(1));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(1));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(0));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(0));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(3));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(3));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(2));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(2));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(13));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(13));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(12));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(12));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(15));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(15));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(14));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(14));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(9));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(9));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(8));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(8));
		returnBuffer.bundle(5).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(11));
		returnBuffer.bundle(5).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(11));
		returnBuffer.bundle(5).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(10));
		returnBuffer.bundle(5).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(10));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(6));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(6));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(7));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(7));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(4));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(4));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(5));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(5));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(2));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(2));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(3));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(3));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(0));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(0));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(1));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(1));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(14));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(14));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(15));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(15));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(12));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(12));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(13));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(13));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(3), stackBuffers[1].bundle(10));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(3), stackBuffers[3].bundle(10));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(11));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(11));
		returnBuffer.bundle(6).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(8));
		returnBuffer.bundle(6).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(8));
		returnBuffer.bundle(6).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(9));
		returnBuffer.bundle(6).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(9));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(11));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(11));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(10));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(10));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(9));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(9));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(8));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(8));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(11), stackBuffers[1].bundle(15));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(11), stackBuffers[3].bundle(15));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(14));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(14));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(13));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(13));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(12));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(12));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(3));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(3));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(2));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(2));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(1));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(1));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(0));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(0));
		returnBuffer.bundle(11).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(7));
		returnBuffer.bundle(11).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(7));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(6));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(6));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(5));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(5));
		returnBuffer.bundle(11).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(4));
		returnBuffer.bundle(11).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(4));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(8));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(8));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(14), stackBuffers[1].bundle(9));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(14), stackBuffers[3].bundle(9));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(10));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(10));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(11));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(11));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(12));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(12));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(13));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(13));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(9), stackBuffers[1].bundle(14));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(9), stackBuffers[3].bundle(14));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(15));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(15));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(0));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(0));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(6), stackBuffers[1].bundle(1));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(6), stackBuffers[3].bundle(1));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(2));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(2));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(3));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(3));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(4));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(4));
		returnBuffer.bundle(8).multAdd(stackBuffers[0].bundle(2), stackBuffers[1].bundle(5));
		returnBuffer.bundle(8).multAdd(stackBuffers[2].bundle(2), stackBuffers[3].bundle(5));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(6));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(6));
		returnBuffer.bundle(8).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(7));
		returnBuffer.bundle(8).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(7));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(9));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(9));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(8));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(8));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(13), stackBuffers[1].bundle(11));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(13), stackBuffers[3].bundle(11));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(12), stackBuffers[1].bundle(10));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(12), stackBuffers[3].bundle(10));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(13));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(13));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(10), stackBuffers[1].bundle(12));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(10), stackBuffers[3].bundle(12));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(15));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(15));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(8), stackBuffers[1].bundle(14));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(8), stackBuffers[3].bundle(14));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(1));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(1));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(0));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(0));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(5), stackBuffers[1].bundle(3));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(5), stackBuffers[3].bundle(3));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(4), stackBuffers[1].bundle(2));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(4), stackBuffers[3].bundle(2));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(5));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(5));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(4));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(4));
		returnBuffer.bundle(9).multSub(stackBuffers[0].bundle(1), stackBuffers[1].bundle(7));
		returnBuffer.bundle(9).multSub(stackBuffers[2].bundle(1), stackBuffers[3].bundle(7));
		returnBuffer.bundle(9).multAdd(stackBuffers[0].bundle(0), stackBuffers[1].bundle(6));
		returnBuffer.bundle(9).multAdd(stackBuffers[2].bundle(0), stackBuffers[3].bundle(6));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(15), stackBuffers[1].bundle(10));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(15), stackBuffers[3].bundle(10));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(14), stackBuffers[1].bundle(11));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(14), stackBuffers[3].bundle(11));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(13), stackBuffers[1].bundle(8));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(13), stackBuffers[3].bundle(8));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(12), stackBuffers[1].bundle(9));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(12), stackBuffers[3].bundle(9));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(11), stackBuffers[1].bundle(14));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(11), stackBuffers[3].bundle(14));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(10), stackBuffers[1].bundle(15));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(10), stackBuffers[3].bundle(15));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(9), stackBuffers[1].bundle(12));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(9), stackBuffers[3].bundle(12));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(8), stackBuffers[1].bundle(13));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(8), stackBuffers[3].bundle(13));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(7), stackBuffers[1].bundle(2));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(7), stackBuffers[3].bundle(2));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(6), stackBuffers[1].bundle(3));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(6), stackBuffers[3].bundle(3));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(5), stackBuffers[1].bundle(0));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(5), stackBuffers[3].bundle(0));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(4), stackBuffers[1].bundle(1));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(4), stackBuffers[3].bundle(1));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(3), stackBuffers[1].bundle(6));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(3), stackBuffers[3].bundle(6));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(2), stackBuffers[1].bundle(7));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(2), stackBuffers[3].bundle(7));
		returnBuffer.bundle(10).multAdd(stackBuffers[0].bundle(1), stackBuffers[1].bundle(4));
		returnBuffer.bundle(10).multAdd(stackBuffers[2].bundle(1), stackBuffers[3].bundle(4));
		returnBuffer.bundle(10).multSub(stackBuffers[0].bundle(0), stackBuffers[1].bundle(5));
		returnBuffer.bundle(10).multSub(stackBuffers[2].bundle(0), stackBuffers[3].bundle(5));
		#pragma endregion
	};

	//propagator bubble
	auto p = [&](const float w1, const float w2) -> float
	{
		return 1.0f / ((w1 + v2->getValue(w1)) *(w2 + v2->getValue(w2)));
	};

	//Katanin contribution, call only once the _flow has been fully calculated and broadcasted
	auto pKataninContribution = [&](const float w1, const float w2) -> float
	{
		float denomW1 = w1 + v2->getValue(w1);
		return static_cast<TRIEffectiveAction *>(_flow)->vertexSingleParticle->getValue(w1) / (denomW1 * denomW1 * (w2 + v2->getValue(w2)));
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
	v4CurrentValue.multAdd(p(cutoff, cutoff + u), buffer1);
	if (u > 2.0f * cutoff)
	{
		integralKernelU(-cutoff, buffer1);
		v4CurrentValue.multAdd(p(cutoff, cutoff - u), buffer1);
	}

	//Katanin contribution
	std::function<void(float, ValueSuperbundle<float, 16> &)> integralKernelSKatanin = [&](float wp, ValueSuperbundle<float, 16> &returnBuffer)->void { integralKernelS(wp, returnBuffer); returnBuffer *= pKataninContribution(wp, s + wp); };
	std::function<void(float, ValueSuperbundle<float, 16> &)> integralKernelTKatanin = [&](float wp, ValueSuperbundle<float, 16> &returnBuffer)->void { integralKernelT(wp, returnBuffer); returnBuffer *= pKataninContribution(wp, t + wp); };
	std::function<void(float, ValueSuperbundle<float, 16> &)> integralKernelUKatanin = [&](float wp, ValueSuperbundle<float, 16> &returnBuffer)->void { integralKernelU(wp, returnBuffer); returnBuffer *= pKataninContribution(wp, u + wp); };

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
	v4CurrentValue /= (2.0f * (float)M_PI);

	for (int b = 0; b < 16; ++b)
	{
		for (int rid = 0; rid < FrgCommon::lattice().size; ++rid) static_cast<TRIEffectiveAction *>(_flow)->vertexTwoParticle->_data[iterator * 16 * FrgCommon::lattice().size + b * FrgCommon::lattice().size + rid] = v4CurrentValue.bundle(b)[rid];
	}
}