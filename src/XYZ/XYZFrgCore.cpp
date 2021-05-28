/**
 * @file XYZFrgCore.cpp
 * @author Finn Lasse Buessen
 * @brief FrgCore implementation for models with diagonal interactions.
 * 
 * @copyright Copyright (c) 2020
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include "lib/InputParser.hpp"
#include "lib/Integrator.hpp"
#include "SpinParser.hpp"
#include "XYZFrgCore.hpp"
#include "XYZEffectiveAction.hpp"

XYZFrgCore::XYZFrgCore(const SpinModel &spinModel, const std::vector<Measurement *> &measurements, const std::map<std::string, std::string> &options) : FrgCore(measurements)
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
	_flowingFunctional = new XYZEffectiveAction(*FrgCommon::cutoff().begin(), spinModel, this);
	_flow = new XYZEffectiveAction();

	//init loadManager
	//stack0
	dataStacks[0] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		&_flowingFunctional->cutoff,
		1);
	//stack1
	dataStacks[1] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexSingleParticle->_data,
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexSingleParticle->size);
	//stack2
	dataStacks[2] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataDD,
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size);
	//stack3
	dataStacks[3] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataXX,
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size);
	//stack4
	dataStacks[4] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataYY,
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size);
	//stack5
	dataStacks[5] = SpinParser::spinParser()->getLoadManager()->addPassiveStack<float>(
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataZZ,
		static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size);
	//stack6
	dataStacks[6] = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		&_flow->cutoff,
		1,
		[&](int) { _flow->cutoff = static_cast<XYZEffectiveAction *>(_flowingFunctional)->cutoff; });
	//stack7
	dataStacks[7] = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		static_cast<XYZEffectiveAction *>(_flow)->vertexSingleParticle->_data,
		static_cast<XYZEffectiveAction *>(_flow)->vertexSingleParticle->size,
		[&](int x) { _calculateVertexSingleParticle(x); },
		1,
		1,
		1);
	//stack8
	dataStacks[8] = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataDD,
		static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->sizeFrequency,
		[&](int x) { _calculateVertexTwoParticle(x); },
		FrgCommon::lattice().size,
		FrgCommon::frequency().size);
	//stack9
	dataStacks[9] = SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataXX,
		static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->sizeFrequency,
		dataStacks[8],
		FrgCommon::lattice().size);
	//stack10
	dataStacks[10] = SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataYY,
		static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->sizeFrequency,
		dataStacks[8],
		FrgCommon::lattice().size);
	//stack11
	dataStacks[11] = SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataZZ,
		static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->sizeFrequency,
		dataStacks[8],
		FrgCommon::lattice().size);
}

XYZFrgCore::~XYZFrgCore()
{
	delete _flowingFunctional;
	delete _flow;
}

void XYZFrgCore::computeStep()
{
	//update cutoff and broadcast
	SpinParser::spinParser()->getLoadManager()->calculate(dataStacks[6]);
	SpinParser::spinParser()->getLoadManager()->broadcast(dataStacks[6]);
	//calculate 1-particle vertices and broadcast (required for Katanin calculation)
	SpinParser::spinParser()->getLoadManager()->calculate(dataStacks[7]);
	SpinParser::spinParser()->getLoadManager()->broadcast(dataStacks[7]);
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
	managedMeasurementStacks.push_back(dataStacks[8]);
	SpinParser::spinParser()->getLoadManager()->calculate(managedMeasurementStacks.data(), int(managedMeasurementStacks.size()));
}

void XYZFrgCore::finalizeStep(float newCutoff)
{
	//determine cutoff set
	float cutoffStep = newCutoff - _flowingFunctional->cutoff;

	//set new cutoff value
	_flowingFunctional->cutoff = newCutoff;

	//add flow to single particle vertex
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexSingleParticle->size; ++i) static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexSingleParticle->_data[i] += cutoffStep * static_cast<XYZEffectiveAction *>(_flow)->vertexSingleParticle->_data[i];

	//add flow to two particle vertex
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size; ++i) static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataDD[i] += cutoffStep * static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataDD[i];
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size; ++i) static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataXX[i] += cutoffStep * static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataXX[i];
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size; ++i) static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataYY[i] += cutoffStep * static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataYY[i];
	#ifndef DISABLE_OMP
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->size; ++i) static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle->_dataZZ[i] += cutoffStep * static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataZZ[i];

	//broadcast updated effective action
	SpinParser::spinParser()->getLoadManager()->broadcast({ dataStacks[0], dataStacks[1], dataStacks[2], dataStacks[3], dataStacks[4], dataStacks[5] });
}

void XYZFrgCore::_calculateVertexSingleParticle(const int iterator)
{
	float cutoff = _flowingFunctional->cutoff;
	XYZVertexSingleParticle *v2 = static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexSingleParticle;
	XYZVertexTwoParticle *v4 = static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle;

	float w;
	v2->expandIterator(iterator, w);
	float v2CurrentValue = 0;

	//term1
	float sum = 0;
	for (auto j = FrgCommon::lattice().getRange(0); j != FrgCommon::lattice().end(); ++j)
	{
		sum += v4->getValue(FrgCommon::lattice().zero(), j, w + cutoff, 0.0f, w - cutoff, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::None);
		sum -= v4->getValue(FrgCommon::lattice().zero(), j, w - cutoff, 0.0f, w + cutoff, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::None);
	}
	v2CurrentValue -= 2.0f * sum;

	//term2
	v2CurrentValue += v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SpinComponent::X, XYZVertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SpinComponent::X, XYZVertexTwoParticle::FrequencyChannel::None);
	v2CurrentValue += v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SpinComponent::Y, XYZVertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SpinComponent::Y, XYZVertexTwoParticle::FrequencyChannel::None);
	v2CurrentValue += v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SpinComponent::Z, XYZVertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SpinComponent::Z, XYZVertexTwoParticle::FrequencyChannel::None);
	v2CurrentValue += v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w + cutoff, w - cutoff, 0.0f, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::None) - v4->getValue(FrgCommon::lattice().zero(), FrgCommon::lattice().zero(), w - cutoff, w + cutoff, 0.0f, SpinComponent::None, XYZVertexTwoParticle::FrequencyChannel::None);

	//prefactor
	v2CurrentValue /= (2.0f * (float)M_PI * (cutoff + v2->getValue(cutoff)));

	static_cast<XYZEffectiveAction *>(_flow)->vertexSingleParticle->_data[iterator] = v2CurrentValue;
}

void XYZFrgCore::_calculateVertexTwoParticle(const int iterator)
{
	float cutoff = _flowingFunctional->cutoff;
	XYZVertexSingleParticle *v2 = static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexSingleParticle;
	XYZVertexTwoParticle *v4 = static_cast<XYZEffectiveAction *>(_flowingFunctional)->vertexTwoParticle;

	float s, t, u;
	v4->expandIterator(iterator, s, t, u);

	//vertex buffers
	ValueSuperbundle<float, 4> buffer1(FrgCommon::lattice().size);
	ValueSuperbundle<float, 4> buffer2(FrgCommon::lattice().size);
	ValueSuperbundle<float, 4> bufferRPA(FrgCommon::lattice().size);
	ValueSuperbundle<float, 4> v4CurrentValue(FrgCommon::lattice().size);
	ValueSuperbundle<float, 4> stackBuffers[4] = {
		ValueSuperbundle<float, 4>(FrgCommon::lattice().size),
		ValueSuperbundle<float, 4>(FrgCommon::lattice().size),
		ValueSuperbundle<float, 4>(FrgCommon::lattice().size),
		ValueSuperbundle<float, 4>(FrgCommon::lattice().size)
	};

	//transfer frequencies
	float w1p = 0.5f * (s + t + u);
	float w1 = 0.5f * (s - t + u);
	float w2p = 0.5f * (s - t - u);
	float w2 = 0.5f * (s + t - u);

	//integrand of the ferquency integral
	auto integralKernelS = [&](const float wp, ValueSuperbundle<float, 4> &returnBuffer) -> void
	{
		//pp-ladder A and B (positive sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab0 = v4->generateAccessBuffer(s, -w1 - wp, -w2 - wp, XYZVertexTwoParticle::FrequencyChannel::S);
		const XYZVertexTwoParticleAccessBuffer<4> ab1 = v4->generateAccessBuffer(s, w1p + wp, -w2p - wp, XYZVertexTwoParticle::FrequencyChannel::S);
		//pp-ladder A and B (positive sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab2 = v4->generateAccessBuffer(s, w2 + wp, w1 + wp, XYZVertexTwoParticle::FrequencyChannel::S);
		const XYZVertexTwoParticleAccessBuffer<4> ab3 = v4->generateAccessBuffer(s, -w2p - wp, w1p + wp, XYZVertexTwoParticle::FrequencyChannel::S);

		v4->getValueSuperbundle(ab0, stackBuffers[0]);
		v4->getValueSuperbundle(ab1, stackBuffers[1]);
		v4->getValueSuperbundle(ab2, stackBuffers[2]);
		v4->getValueSuperbundle(ab3, stackBuffers[3]);

		//calculate flow
		returnBuffer.reset();

		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));

		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));

		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));

		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));
	};

	auto integralKernelT = [&](const float wp, ValueSuperbundle<float, 4> &returnBuffer) -> void
	{
		//RPA diagram A and B equal chalice diagram A and inverse chalice diagram B, respectively (negative sign)
		//chalice diagram A (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab0 = v4->generateAccessBuffer(w1 - wp, t, w1p + wp, XYZVertexTwoParticle::FrequencyChannel::T);
		//inverse chalice diagram B (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab1 = v4->generateAccessBuffer(w2p - wp, t, -w2 - wp, XYZVertexTwoParticle::FrequencyChannel::T);
		//chalice diagram A (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab2 = v4->generateAccessBuffer(w1p + wp, t, w1 - wp, XYZVertexTwoParticle::FrequencyChannel::T);
		//inverse chalice diagram B (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab3 = v4->generateAccessBuffer(w2 + wp, t, wp - w2p, XYZVertexTwoParticle::FrequencyChannel::T);

		v4->getValueSuperbundle(ab0, stackBuffers[0]);
		v4->getValueSuperbundle(ab1, stackBuffers[1]);
		v4->getValueSuperbundle(ab2, stackBuffers[2]);
		v4->getValueSuperbundle(ab3, stackBuffers[3]);

		//calculate flow
		returnBuffer.reset();
		bufferRPA.reset();
		for (int rid = 0; rid < FrgCommon::lattice().size; ++rid)
		{
			//lattice bubble
			const LatticeOverlap &overlap = FrgCommon::lattice().getOverlap(rid);

			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(static_cast<int>(SpinComponent::X))[rid] += stackBuffers[0].bundle(static_cast<int>(overlap.transformedX1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(static_cast<int>(overlap.transformedX2[i]))[overlap.rid2[i]];
			}

			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(static_cast<int>(SpinComponent::Y))[rid] += stackBuffers[0].bundle(static_cast<int>(overlap.transformedY1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(static_cast<int>(overlap.transformedY2[i]))[overlap.rid2[i]];
			}

			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(static_cast<int>(SpinComponent::Z))[rid] += stackBuffers[0].bundle(static_cast<int>(overlap.transformedZ1[i]))[overlap.rid1[i]] * stackBuffers[1].bundle(static_cast<int>(overlap.transformedZ2[i]))[overlap.rid2[i]];
			}

			for (int i = 0; i < overlap.size; ++i)
			{
				bufferRPA.bundle(static_cast<int>(SpinComponent::None))[rid] += stackBuffers[0].bundle(3)[overlap.rid1[i]] * stackBuffers[1].bundle(3)[overlap.rid2[i]];
			}
		}
		returnBuffer.multAdd(4.0f, bufferRPA);

		//chalice diagram B (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab4 = v4->generateAccessBuffer(w2p - wp, -w2 - wp, t, XYZVertexTwoParticle::FrequencyChannel::U);
		//inverse chalice diagram A (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab5 = v4->generateAccessBuffer(w1 - wp, -w1p - wp, -t, XYZVertexTwoParticle::FrequencyChannel::U);
		//chalice diagram B (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab6 = v4->generateAccessBuffer(w2 + wp, wp - w2p, t, XYZVertexTwoParticle::FrequencyChannel::U);
		//inverse chalice diagram A (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab7 = v4->generateAccessBuffer(w1p + wp, wp - w1, -t, XYZVertexTwoParticle::FrequencyChannel::U);

		const float valCbx = v4->getValueLocal(SpinComponent::X, ab4);
		const float valCby = v4->getValueLocal(SpinComponent::Y, ab4);
		const float valCbz = v4->getValueLocal(SpinComponent::Z, ab4);
		const float valCbd = v4->getValueLocal(SpinComponent::None, ab4);
		const float valICax = v4->getValueLocal(SpinComponent::X, ab5);
		const float valICay = v4->getValueLocal(SpinComponent::Y, ab5);
		const float valICaz = v4->getValueLocal(SpinComponent::Z, ab5);
		const float valICad = v4->getValueLocal(SpinComponent::None, ab5);
		const float valCbx2 = v4->getValueLocal(SpinComponent::X, ab6);
		const float valCby2 = v4->getValueLocal(SpinComponent::Y, ab6);
		const float valCbz2 = v4->getValueLocal(SpinComponent::Z, ab6);
		const float valCbd2 = v4->getValueLocal(SpinComponent::None, ab6);
		const float valICax2 = v4->getValueLocal(SpinComponent::X, ab7);
		const float valICay2 = v4->getValueLocal(SpinComponent::Y, ab7);
		const float valICaz2 = v4->getValueLocal(SpinComponent::Z, ab7);
		const float valICad2 = v4->getValueLocal(SpinComponent::None, ab7);

		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), valCbd);
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), valCby);
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), valCbz);
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), valCbx);
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(valICad, stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(valICaz, stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(valICay, stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(valICax, stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), valCbd2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), valCby2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), valCbz2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), valCbx2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(valICad2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(valICaz2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multAdd(valICay2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(valICax2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));

		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), valCbd);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), valCbx);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), valCbz);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), valCby);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(valICad, stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(valICaz, stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(valICax, stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(valICay, stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), valCbd2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), valCbx2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), valCbz2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), valCby2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(valICad2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(valICaz2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multAdd(valICax2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(valICay2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));

		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), valCbd);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), valCbx);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), valCby);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), valCbz);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(valICad, stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(valICax, stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(valICay, stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(valICaz, stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), valCbd2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), valCbx2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), valCby2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), valCbz2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(valICad2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(valICax2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multAdd(valICay2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(valICaz2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));

		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), valCbd);
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), valCbx);
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), valCby);
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), valCbz);
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(valICad, stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(valICax, stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(valICay, stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(valICaz, stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), valCbd2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), valCbx2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), valCby2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), valCbz2);
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(valICad2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(valICax2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(valICay2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(valICaz2, stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
	};

	auto integralKernelU = [&](const float wp, ValueSuperbundle<float, 4> &returnBuffer) -> void
	{
		//u-Channel, to be combined with P(wp, u + wp) + P(u + wp, wp)
		//ph-ladder A and B, respectively (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab0 = v4->generateAccessBuffer(w1 + wp, wp - w2p, u, XYZVertexTwoParticle::FrequencyChannel::U);
		const XYZVertexTwoParticleAccessBuffer<4> ab1 = v4->generateAccessBuffer(w1p + wp, w2 - wp, u, XYZVertexTwoParticle::FrequencyChannel::U);
		//ph-ladder A and B, respectively (negative sign)
		const XYZVertexTwoParticleAccessBuffer<4> ab2 = v4->generateAccessBuffer(w2p - wp, -w1 - wp, u, XYZVertexTwoParticle::FrequencyChannel::U);
		const XYZVertexTwoParticleAccessBuffer<4> ab3 = v4->generateAccessBuffer(w2 - wp, w1p + wp, u, XYZVertexTwoParticle::FrequencyChannel::U);

		v4->getValueSuperbundle(ab0, stackBuffers[0]);
		v4->getValueSuperbundle(ab1, stackBuffers[1]);
		v4->getValueSuperbundle(ab2, stackBuffers[2]);
		v4->getValueSuperbundle(ab3, stackBuffers[3]);

		//calculate flow
		returnBuffer.reset();

		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::X)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));

		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Y)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));

		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::Z)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));

		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[0].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[1].bundle(static_cast<int>(SpinComponent::Z)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::None)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::None)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::X)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::X)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Y)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Y)));
		returnBuffer.bundle(static_cast<int>(SpinComponent::None)).multSub(stackBuffers[2].bundle(static_cast<int>(SpinComponent::Z)), stackBuffers[3].bundle(static_cast<int>(SpinComponent::Z)));
	};

	//propagator bubble
	auto p = [&](const float w1, const float w2) -> float
	{
		return 1.0f / ((w1 + v2->getValue(w1)) *(w2 + v2->getValue(w2)));
	};

	//Katanin contribution, call only once the flow has been fully calculated and broadcasted
	auto pKataninContribution = [&](const float w1, const float w2) -> float
	{
		float denomW1 = w1 + v2->getValue(w1);
		return static_cast<XYZEffectiveAction *>(_flow)->vertexSingleParticle->getValue(w1) / (denomW1 * denomW1 * (w2 + v2->getValue(w2)));
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
	std::function<void(float, ValueSuperbundle<float, 4> &)> integralKernelSKatanin = [&](float wp, ValueSuperbundle<float, 4> &returnBuffer)->void { integralKernelS(wp, returnBuffer); returnBuffer *= pKataninContribution(wp, s + wp); };
	std::function<void(float, ValueSuperbundle<float, 4> &)> integralKernelTKatanin = [&](float wp, ValueSuperbundle<float, 4> &returnBuffer)->void { integralKernelT(wp, returnBuffer); returnBuffer *= pKataninContribution(wp, t + wp); };
	std::function<void(float, ValueSuperbundle<float, 4> &)> integralKernelUKatanin = [&](float wp, ValueSuperbundle<float, 4> &returnBuffer)->void { integralKernelU(wp, returnBuffer); returnBuffer *= pKataninContribution(wp, u + wp); };

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

	for (int rid = 0; rid < FrgCommon::lattice().size; ++rid) static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataXX[iterator * FrgCommon::lattice().size + rid] = v4CurrentValue.bundle(static_cast<int>(SpinComponent::X))[rid];
	for (int rid = 0; rid < FrgCommon::lattice().size; ++rid) static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataYY[iterator * FrgCommon::lattice().size + rid] = v4CurrentValue.bundle(static_cast<int>(SpinComponent::Y))[rid];
	for (int rid = 0; rid < FrgCommon::lattice().size; ++rid) static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataZZ[iterator * FrgCommon::lattice().size + rid] = v4CurrentValue.bundle(static_cast<int>(SpinComponent::Z))[rid];
	for (int rid = 0; rid < FrgCommon::lattice().size; ++rid) static_cast<XYZEffectiveAction *>(_flow)->vertexTwoParticle->_dataDD[iterator * FrgCommon::lattice().size + rid] = v4CurrentValue.bundle(static_cast<int>(SpinComponent::None))[rid];
}