/**
 * @file TRIVertexTwoParticle.hpp
 * @author Finn Lasse Buessen
 * @brief Two-particle vertex implementation for time reversal invariant models.
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <istream>
#include "lib/ValueBundle.hpp"
#include "lib/Assert.hpp"
#include "FrgCommon.hpp"

/**
 * @brief Buffer for frequency interpolation information. 
 * @details The buffer contains a list of memory offsets (number of elements) in the frequency dimensions of the two-particle vertex that correspond to all the support values that are part of the interpolation. 
 * The memory offset in frequency space must be complemented with the linear memory offset in the lattice site dimension. 
 * Each support value is assigned a weight factor. 
 * Each support value is also complemented by a sign factor that is either +1 or -1 (and depends on the two vertex channel arguments) and by the information whether a lattice site exchange should be performed upon accessing the support value. 
 * 
 * @tparam size Number of support sites for the interpolation. 
 */
template <int size> struct TRIVertexTwoParticleAccessBuffer
{
public:
	/**
	 * @brief Construct a new TRIVertexTwoParticleAccessBuffer object and initialize the sign factors to plus one. 
	 */
	TRIVertexTwoParticleAccessBuffer() : pairExchange(false) 
	{
		for (int i = 0; i < 16 * size; ++i) (&sign[0][0][0])[i] = 1.0f;
	}

	int frequencyOffsets[size]; ///< Linear memory offset (number of elements) in the frequency dimensions of the two-particle vertex.  
	float frequencyWeights[size]; ///< Weight factors of the support values. 
	float sign[size][4][4]; ///< Sign factors of the support values. 
	bool pairExchange; ///< Site exchange indicator. 
};

/**
 * @brief Two-particle vertex implementation for time reversal invariant models.
 */
struct TRIVertexTwoParticle
{
public:
	/**
	 * @brief Indicator for frequency channels that exactly lie on frequency mesh points.  
	 */
	enum struct FrequencyChannel
	{
		S, ///< s-channel. 
		T, ///< t-channel. 
		U, ///< u-channel. 
		All, ///< All channels. 
		None ///< No channel. 
	};

	/**
	 * @brief Construct a new TRIVertexTwoParticle object and initialize all entries to zero. 
	 */
	TRIVertexTwoParticle()
	{
		//store width in all memory dimensions
		_memoryStep[3] = FrgCommon::lattice().size;
		_memoryStep[2] = 4 * _memoryStep[3];
		_memoryStep[1] = 4 * _memoryStep[2];
		_memoryStep[0] = FrgCommon::frequency().size * _memoryStep[1];

		sizeFrequency = FrgCommon::frequency().size * FrgCommon::frequency().size * (FrgCommon::frequency().size + 1) / 2;
		size = 16 * FrgCommon::lattice().size * sizeFrequency;

		//alloc and init memory
		_data = new float[size];
		memset(_data, 0, sizeof(float) * size);
	}

	/**
	 * @brief Destroy the TRIVertexTwoParticle object. 
	 */
	~TRIVertexTwoParticle()
	{
		delete[] _data;
	}

	/**
	 * @brief Expand a linear iterator in the range [0,size) that iterates over all vertex entries. 
	 * 
	 * @param[in] iterator Linear iterator. 
	 * @param[out] i1 Lattice site iterator. 
	 * @param[out] s First frequency argument. 
	 * @param[out] t Second frequency argument. 
	 * @param[out] u Third frequency argument. 
	 * @param[out] s1 First vertex channel. 
	 * @param[out] s2 Second vertex channel. 
	 */
	void expandIterator(int iterator, LatticeIterator &i1, float &s, float &t, float &u, SpinComponent &s1, SpinComponent &s2) const
	{
		ASSERT(iterator >= 0 && iterator < size);
		ASSERT(&s != &t);
		ASSERT(&t != &u);
		ASSERT(&s != &u);
		ASSERT(&s1 != &s2);

		int it = iterator;

		int su = it / _memoryStep[0];
		it = it % _memoryStep[0];
		t = FrgCommon::frequency()._data[it / _memoryStep[1]];
		it = it % _memoryStep[1];
		s1 = static_cast<SpinComponent>(it / _memoryStep[2]);
		it = it % _memoryStep[2];
		s2 = static_cast<SpinComponent>(it / _memoryStep[3]);
		it = it % _memoryStep[3];
		i1 = FrgCommon::lattice().fromParametrization(it);

		for (int so = 0; so <= su; ++so)
		{
			for (int uo = 0; uo <= so; ++uo)
			{
				if (su == so * (so + 1) / 2 + uo)
				{
					s = FrgCommon::frequency()._data[so];
					u = FrgCommon::frequency()._data[uo];
					return;
				}
			}
		}
	}

	/**
	 * @brief Expand a linear iterator in the range [0,sizeFrequency) that iterates over all paramtetrized frequency values. 
	 * 
	 * @param[in] iterator Linear iterator. 
	 * @param[out] s First frequency argument. 
	 * @param[out] t Second frequency argument. 
	 * @param[out] u Third frequency argument. 
	 */
	void expandIterator(int iterator, float &s, float &t, float &u) const
	{
		ASSERT(iterator >= 0 && iterator < sizeFrequency);
		ASSERT(iterator >= 0 && iterator < size);
		ASSERT(&s != &t);
		ASSERT(&t != &u);
		ASSERT(&s != &u);

		int su = iterator / FrgCommon::frequency().size;
		t = FrgCommon::frequency()._data[iterator % FrgCommon::frequency().size];

		for (int so = 0; so <= su; ++so)
		{
			for (int uo = 0; uo <= so; ++uo)
			{
				if (su == so * (so + 1) / 2 + uo)
				{
					s = FrgCommon::frequency()._data[so];
					u = FrgCommon::frequency()._data[uo];
					return;
				}
			}
		}
	}

	/**
	 * @brief Directly access a vertex value via a linear iterator in the range [0,size). 
	 * 
	 * @param iterator Linear iterator. 
	 * @return float& Vertex value. 
	 */
	float &getValueRef(const int iterator) const
	{
		ASSERT(iterator >= 0 && iterator < size);

		return _data[iterator];
	}

	/**
	 * @brief Access vertex value at arbitrary lattice sites, frequencies, and spin components. 
	 * 
	 * @param i1 First lattice site argument. 
	 * @param i2 Second lattice site argument. 
	 * @param s First frequency argument. 
	 * @param t Second frequency argument. 
	 * @param u Third frequency argument. 
	 * @param s1 First vertex channel. 
	 * @param s2 Second vertex channel. 
	 * @param channel Frequency channel. 
	 * @return float Vertex value. 
	 */
	float getValue(LatticeIterator i1, LatticeIterator i2, float s, float t, float u, SpinComponent s1, SpinComponent s2, const FrequencyChannel channel) const
	{
		ASSERT(channel == FrequencyChannel::S || channel == FrequencyChannel::T || channel == FrequencyChannel::U || channel == FrequencyChannel::None || channel == FrequencyChannel::All);

		//map to positive frequency sector
		float sign = 1.0f;
		if (s < 0)
		{
			s = -s;
			std::swap(s1, s2);
			std::swap(i1, i2);
		}
		if (t < 0)
		{
			t = -t;
			sign *= _zeta(static_cast<int>(s1)) *_zeta(static_cast<int>(s2));
		}
		if (u < 0)
		{
			u = -u;
			std::swap(s1, s2);
			std::swap(i1, i2);
			sign *= _zeta(static_cast<int>(s1)) *_zeta(static_cast<int>(s2));
		}

		int siteOffset = FrgCommon::lattice().symmetryTransform(i1, i2, s1, s2);

		if (channel == FrequencyChannel::S)
		{
			int exactS = FrgCommon::frequency().offset(s);

			int lowerT, upperT;
			float biasT;
			int lowerU, upperU;
			float biasU;

			FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);
			FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

			return sign * (
				(1 - biasU) * (
				(1 - biasT) * (_directAccessMapFrequencyExchange(siteOffset, exactS, lowerT, lowerU, s1, s2)) + biasT * (_directAccessMapFrequencyExchange(siteOffset, exactS, upperT, lowerU, s1, s2))
					) + biasU * (
					(1 - biasT) * (_directAccessMapFrequencyExchange(siteOffset, exactS, lowerT, upperU, s1, s2)) + biasT * (_directAccessMapFrequencyExchange(siteOffset, exactS, upperT, upperU, s1, s2))
						));
		}
		else if (channel == FrequencyChannel::T)
		{
			int exactT = FrgCommon::frequency().offset(t);

			int lowerS, upperS;
			float biasS;
			int lowerU, upperU;
			float biasU;

			FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
			FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

			return sign * (
				(1 - biasU) * (
				(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, exactT, lowerU, s1, s2)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, exactT, lowerU, s1, s2))
					) + biasU * (
					(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, exactT, upperU, s1, s2)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, exactT, upperU, s1, s2))
						));
		}
		else if (channel == FrequencyChannel::U)
		{
			int exactU = FrgCommon::frequency().offset(u);

			int lowerS, upperS;
			float biasS;
			int lowerT, upperT;
			float biasT;

			FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
			FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);

			return sign * (
				(1 - biasT) * (
				(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, lowerT, exactU, s1, s2)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, lowerT, exactU, s1, s2))
					) + biasT * (
					(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, upperT, exactU, s1, s2)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, upperT, exactU, s1, s2))
						));
		}
		else if (channel == FrequencyChannel::None)
		{
			int lowerS, upperS;
			float biasS;
			int lowerT, upperT;
			float biasT;
			int lowerU, upperU;
			float biasU;

			FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
			FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);
			FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

			return sign * (
				(1 - biasU) * (
				(1 - biasT) * (
					(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, lowerT, lowerU, s1, s2)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, lowerT, lowerU, s1, s2))
					) + biasT * (
					(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, upperT, lowerU, s1, s2)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, upperT, lowerU, s1, s2))
						)
					) + biasU * (
					(1 - biasT) * (
						(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, lowerT, upperU, s1, s2)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, lowerT, upperU, s1, s2))
						) + biasT * (
						(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, upperT, upperU, s1, s2)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, upperT, upperU, s1, s2))
							)
						));
		}
		else if (channel == FrequencyChannel::All)
		{
			int exactS = FrgCommon::frequency().offset(s);
			int exactT = FrgCommon::frequency().offset(t);
			int exactU = FrgCommon::frequency().offset(u);
			return sign * _directAccessMapFrequencyExchange(siteOffset, exactS, exactT, exactU, s1, s2);
		}
		else throw Exception(Exception::Type::ArgumentError, "Specified frequency channel does not exist");
	}

	/**
	 * @brief Access vertex value at arbitrary lattice sites and spin components via a given access buffer.
	 * 
	 * @tparam n Number of support sites in the access buffer. 
	 * @param i1 First lattice site argument. 
	 * @param i2 Second lattice site argument. 
	 * @param s1 First vertex channel. 
	 * @param s2 Second vertex channel. 
	 * @param accessBuffer Access buffer. 
	 * @return float Vertex value. 
	 */
	template <int n> float getValue(LatticeIterator i1, LatticeIterator i2, SpinComponent s1, SpinComponent s2, const TRIVertexTwoParticleAccessBuffer<n> &accessBuffer) const
	{
		if (accessBuffer.pairExchange)
		{
			std::swap(i1, i2);
			std::swap(s1, s2);
		}
		int siteOffset = FrgCommon::lattice().symmetryTransform(i1, i2, s1, s2);
		int spinOffset = (4 * static_cast<int>(s1) + static_cast<int>(s2)) *FrgCommon::lattice().size;

		float value = 0.0f;
		for (int i = 0; i < n; ++i) value += accessBuffer.sign[i][static_cast<int>(s1)][static_cast<int>(s2)] * accessBuffer.frequencyWeights[i] * _data[accessBuffer.frequencyOffsets[i] + spinOffset + siteOffset];
		return value;
	}

	/**
	 * @brief Access vertex value locally at arbitrary spin components via a given access buffer. 
	 * 
	 * @tparam n Number of support sites in the access buffer. 
	 * @param s1 First vertex channel. 
	 * @param s2 Second vertex channel. 
	 * @param accessBuffer Access buffer. 
	 * @return float Vertex value. 
	 */
	template <int n> float getValueLocal(SpinComponent s1, SpinComponent s2, const TRIVertexTwoParticleAccessBuffer<n> &accessBuffer) const
	{
		if (accessBuffer.pairExchange) std::swap(s1, s2);
		int spinOffset = (4 * static_cast<int>(s1) + static_cast<int>(s2)) *FrgCommon::lattice().size;

		float value = 0.0f;
		for (int i = 0; i < n; ++i) value += accessBuffer.sign[i][static_cast<int>(s1)][static_cast<int>(s2)] * accessBuffer.frequencyWeights[i] * _data[accessBuffer.frequencyOffsets[i] + spinOffset];
		return value;
	}

	/**
	 * @brief Bundled vertex access on all lattice sites and spin components simultaneously via a given access buffer. 
	 * 
	 * @tparam n Number of support sites in the access buffer. 
	 * @param[in] accessBuffer Access buffer. 
	 * @param[out] superbundle Vertex value bundle. 
	 */
	template <int n> void getValueSuperbundle(const TRIVertexTwoParticleAccessBuffer<n> &accessBuffer, ValueSuperbundle<float, 16> &superbundle) const
	{
		ASSERT(superbundle.bundle(0).size() == FrgCommon::lattice().size);

		superbundle.reset();
		const LatticeSiteDescriptor *sites = (accessBuffer.pairExchange) ? FrgCommon::lattice().getInvertedSites() : FrgCommon::lattice().getSites();

		for (int i = 0; i < n; ++i)
		{
			for (int s1 = 0; s1 < 4; ++s1)
			{
				for (int s2 = 0; s2 < 4; ++s2)
				{
					for (int j = 0; j < FrgCommon::lattice().size; ++j)
					{
						int s1t = (accessBuffer.pairExchange) ? s2 : s1;
						int s2t = (accessBuffer.pairExchange) ? s1 : s2;
						if (s1t < 3) s1t = static_cast<int>(sites[j].spinPermutation[s1t]);
						if (s2t < 3) s2t = static_cast<int>(sites[j].spinPermutation[s2t]);
						int spinOffset = (4 * s1t + s2t) * FrgCommon::lattice().size;

						superbundle.bundle(4 * s1 + s2)[j] += accessBuffer.sign[i][s1][s2] * accessBuffer.frequencyWeights[i] * _data[accessBuffer.frequencyOffsets[i] + spinOffset + sites[j].rid];
					}
				}
			}
		}
	}

	/**
	 * @brief Generate an access buffer for a set of frequencies where one of them (specified by channel) exactly lies on the frequency mesh. 
	 * 
	 * @param s First frequency argument. 
	 * @param t Second frequency argument. 
	 * @param u Third frequency argument. 
	 * @param channel Frequency channel. 
	 * @return TRIVertexTwoParticleAccessBuffer<4> Access buffer. 
	 */
	TRIVertexTwoParticleAccessBuffer<4> generateAccessBuffer(float s, float t, float u, const FrequencyChannel channel) const
	{
		ASSERT(channel == FrequencyChannel::S || channel == FrequencyChannel::T || channel == FrequencyChannel::U);

		TRIVertexTwoParticleAccessBuffer<4> ab;

		//map to positive frequency sector
		if (s < 0)
		{
			s = -s;
			ab.pairExchange = !ab.pairExchange;
		}
		if (t < 0)
		{
			t = -t;
			for (int i = 0; i < 4; ++i)
			{
				for (int s1 = 0; s1 <= 3; ++s1)
				{
					for (int s2 = 0; s2 <= 3; ++s2)
					{
						ab.sign[i][s1][s2] *= _zeta(s1) * _zeta(s2);
					}
				}
			}
		}
		if (u < 0)
		{
			u = -u;
			ab.pairExchange = !ab.pairExchange;
			for (int i = 0; i < 4; ++i)
			{
				for (int s1 = 0; s1 <= 3; ++s1)
				{
					for (int s2 = 0; s2 <= 3; ++s2)
					{
						ab.sign[i][s1][s2] *= _zeta(s1) * _zeta(s2);
					}
				}
			}
		}

		//interpolate frequency
		if (channel == FrequencyChannel::S)
		{
			int exactS = FrgCommon::frequency().offset(s);
			int lowerT, upperT;
			float biasT;
			int lowerU, upperU;
			float biasU;
			FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);
			FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

			ab.frequencyWeights[0] = (1 - biasU) * (1 - biasT);
			_generateAccessBufferOffsetMapFrequencyExchange(exactS, lowerT, lowerU, 0, ab);
			ab.frequencyWeights[1] = (1 - biasU) * biasT;
			_generateAccessBufferOffsetMapFrequencyExchange(exactS, upperT, lowerU, 1, ab);
			ab.frequencyWeights[2] = biasU * (1 - biasT);
			_generateAccessBufferOffsetMapFrequencyExchange(exactS, lowerT, upperU, 2, ab);
			ab.frequencyWeights[3] = biasU * biasT;
			_generateAccessBufferOffsetMapFrequencyExchange(exactS, upperT, upperU, 3, ab);
		}
		else if (channel == FrequencyChannel::T)
		{
			int exactT = FrgCommon::frequency().offset(t);
			int lowerS, upperS;
			float biasS;
			int lowerU, upperU;
			float biasU;
			FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
			FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

			ab.frequencyWeights[0] = (1 - biasU) * (1 - biasS);
			_generateAccessBufferOffsetMapFrequencyExchange(lowerS, exactT, lowerU, 0, ab);
			ab.frequencyWeights[1] = (1 - biasU) * biasS;
			_generateAccessBufferOffsetMapFrequencyExchange(upperS, exactT, lowerU, 1, ab);
			ab.frequencyWeights[2] = biasU * (1 - biasS);
			_generateAccessBufferOffsetMapFrequencyExchange(lowerS, exactT, upperU, 2, ab);
			ab.frequencyWeights[3] = biasU * biasS;
			_generateAccessBufferOffsetMapFrequencyExchange(upperS, exactT, upperU, 3, ab);

		}
		else if (channel == FrequencyChannel::U)
		{
			int exactU = FrgCommon::frequency().offset(u);
			int lowerS, upperS;
			float biasS;
			int lowerT, upperT;
			float biasT;
			FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
			FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);

			ab.frequencyWeights[0] = (1 - biasT) * (1 - biasS);
			_generateAccessBufferOffsetMapFrequencyExchange(lowerS, lowerT, exactU, 0, ab);
			ab.frequencyWeights[1] = (1 - biasT) * biasS;
			_generateAccessBufferOffsetMapFrequencyExchange(upperS, lowerT, exactU, 1, ab);
			ab.frequencyWeights[2] = biasT * (1 - biasS);
			_generateAccessBufferOffsetMapFrequencyExchange(lowerS, upperT, exactU, 2, ab);
			ab.frequencyWeights[3] = biasT * biasS;
			_generateAccessBufferOffsetMapFrequencyExchange(upperS, upperT, exactU, 3, ab);
		}

		return ab;
	}

	/**
	 * @brief Generate an access buffer for an arbitrary set of frequencies. 
	 * 
	 * @param s First frequency argument. 
	 * @param t Second frequency argument. 
	 * @param u Third frequency argument. 
	 * @return TRIVertexTwoParticleAccessBuffer<8> Access buffer. 
	 */
	TRIVertexTwoParticleAccessBuffer<8> generateAccessBuffer(float s, float t, float u) const
	{
		TRIVertexTwoParticleAccessBuffer<8> ab;

		//map to positive frequency sector
		if (s < 0)
		{
			s = -s;
			ab.pairExchange = !ab.pairExchange;
		}
		if (t < 0)
		{
			t = -t;
			for (int i = 0; i < 8; ++i)
			{
				for (int s1 = 0; s1 <= 3; ++s1)
				{
					for (int s2 = 0; s2 <= 3; ++s2)
					{
						if (ab.pairExchange) ab.sign[i][s1][s2] *= _zeta(s2) * _zeta(s1);
						else ab.sign[i][s1][s2] *= _zeta(s1) * _zeta(s2);
					}
				}
			}
		}
		if (u < 0)
		{
			u = -u;
			ab.pairExchange = !ab.pairExchange;
			for (int i = 0; i < 8; ++i)
			{
				for (int s1 = 0; s1 <= 3; ++s1)
				{
					for (int s2 = 0; s2 <= 3; ++s2)
					{
						if (ab.pairExchange) ab.sign[i][s1][s2] *= _zeta(s2) * _zeta(s1);
						else ab.sign[i][s1][s2] *= _zeta(s1) * _zeta(s2);
					}
				}
			}
		}

		int lowerS, upperS;
		float biasS;
		int lowerT, upperT;
		float biasT;
		int lowerU, upperU;
		float biasU;

		FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
		FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);
		FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

		ab.frequencyWeights[0] = (1 - biasT) * (1 - biasS) * (1 - biasU);
		_generateAccessBufferOffsetMapFrequencyExchange(lowerS, lowerT, lowerU, 0, ab);
		ab.frequencyWeights[1] = (1 - biasT) * biasS * (1 - biasU);
		_generateAccessBufferOffsetMapFrequencyExchange(upperS, lowerT, lowerU, 1, ab);
		ab.frequencyWeights[2] = biasT * (1 - biasS) * (1 - biasU);
		_generateAccessBufferOffsetMapFrequencyExchange(lowerS, upperT, lowerU, 2, ab);
		ab.frequencyWeights[3] = biasT * biasS * (1 - biasU);
		_generateAccessBufferOffsetMapFrequencyExchange(upperS, upperT, lowerU, 3, ab);
		ab.frequencyWeights[4] = (1 - biasT) * (1 - biasS) * biasU;
		_generateAccessBufferOffsetMapFrequencyExchange(lowerS, lowerT, upperU, 4, ab);
		ab.frequencyWeights[5] = (1 - biasT) * biasS * biasU;
		_generateAccessBufferOffsetMapFrequencyExchange(upperS, lowerT, upperU, 5, ab);
		ab.frequencyWeights[6] = biasT * (1 - biasS) * biasU;
		_generateAccessBufferOffsetMapFrequencyExchange(lowerS, upperT, upperU, 6, ab);
		ab.frequencyWeights[7] = biasT * biasS * biasU;
		_generateAccessBufferOffsetMapFrequencyExchange(upperS, upperT, upperU, 7, ab);

		return ab;
	}

	/**
	 * @brief Directly access a vertex via given spin components, frequency, and site offsets, where sOffset may be lesser than uOffset.
	 * 
	 * @param siteOffset Lattice site offset (number of elements). 
	 * @param sOffset First frequency offset (number of elements). 
	 * @param tOffset Second frequency offset (number of elements). 
	 * @param uOffset Third frequency offset (number of elements). 
	 * @param s1 First vertex channel.
	 * @param s2 Second vertex channel. 
	 * @return float Vertex value. 
	 */
	float _directAccessMapFrequencyExchange(const int siteOffset, const int sOffset, const int tOffset, const int uOffset, const SpinComponent s1, const SpinComponent s2) const
	{
		ASSERT(siteOffset >= 0 && siteOffset < FrgCommon::lattice().size);
		ASSERT(sOffset >= 0 && sOffset < FrgCommon::frequency().size);
		ASSERT(tOffset >= 0 && tOffset < FrgCommon::frequency().size);
		ASSERT(uOffset >= 0 && uOffset < FrgCommon::frequency().size);

		if (sOffset < uOffset) return -_zeta(static_cast<int>(s2)) *_directAccess(siteOffset, uOffset, tOffset, sOffset, s1, s2);
		else return _directAccess(siteOffset, sOffset, tOffset, uOffset, s1, s2);
	}

	/**
	 * @brief Directly access a vertex via given spin components, frequency, and site offsets, where sOffset >= uOffset.
	 * 
	 * @param siteOffset Lattice site offset (number of elements). 
	 * @param sOffset First frequency offset (number of elements). 
	 * @param tOffset Second frequency offset (number of elements). 
	 * @param uOffset Third frequency offset (number of elements). 
	 * @param s1 First vertex channel.
	 * @param s2 Second vertex channel. 
	 * @return float Vertex value. 
	 */
	float _directAccess(const int siteOffset, const int sOffset, const int tOffset, const int uOffset, const SpinComponent s1, const SpinComponent s2) const
	{
		ASSERT(siteOffset >= 0 && siteOffset < FrgCommon::lattice().size);
		ASSERT(sOffset >= 0 && sOffset < FrgCommon::frequency().size);
		ASSERT(tOffset >= 0 && tOffset < FrgCommon::frequency().size);
		ASSERT(uOffset >= 0 && uOffset < FrgCommon::frequency().size);
		ASSERT(sOffset >= uOffset);

		int suOffset = sOffset * (sOffset + 1) / 2 + uOffset;
		return _data[suOffset * _memoryStep[0] + tOffset * _memoryStep[1] + static_cast<int>(s1) *_memoryStep[2] + static_cast<int>(s2) *_memoryStep[3] + siteOffset];
	}

	/**
	 * @brief Calculate the total memory offset from given frequency offsets, where sOffset may be lesser than uOffset. 
	 * Write the result to the specified entry in an access buffer. 
	 * 
	 * @tparam n Number of support sites in the access buffer. 
	 * @param[in] sOffset First frequency offset (number of elements). 
	 * @param[in] tOffset Second frequency offset (number of elements). 
	 * @param[in] uOffset Third frequency offset (number of elements). 
	 * @param[in] abIndex Index of the access buffer entry to write. 
	 * @param ab Access buffer. 
	 */
	template <int n> void _generateAccessBufferOffsetMapFrequencyExchange(const int sOffset, const int tOffset, const int uOffset, const int abIndex, TRIVertexTwoParticleAccessBuffer<n> &ab) const
	{
		ASSERT(sOffset >= 0 && sOffset < FrgCommon::frequency().size);
		ASSERT(tOffset >= 0 && tOffset < FrgCommon::frequency().size);
		ASSERT(uOffset >= 0 && uOffset < FrgCommon::frequency().size);
		ASSERT(abIndex >= 0 && abIndex < n);

		if (sOffset < uOffset)
		{
			for (int s1 = 0; s1 <= 3; ++s1)
			{
				for (int s2 = 0; s2 <= 3; ++s2)
				{
					if (ab.pairExchange) ab.sign[abIndex][s1][s2] *= -_zeta(s1);
					else ab.sign[abIndex][s1][s2] *= -_zeta(s2);
				}
			}
			ab.frequencyOffsets[abIndex] = (uOffset * (uOffset + 1) / 2 + sOffset) * _memoryStep[0] + tOffset * _memoryStep[1];
		}
		else
		{
			ab.frequencyOffsets[abIndex] = (sOffset * (sOffset + 1) / 2 + uOffset) * _memoryStep[0] + tOffset * _memoryStep[1];
		}
	}

	/**
	 * @brief Helper function for vertex symmetries. Returns -1 if the argument is spin-like and +1 if it is density-like. 
	 * 
	 * @param s1 Spin argument. 
	 * @return float Function value. 
	 */
	float _zeta(const int s1) const
	{
		return (s1 <= 2) ? -1.0f : 1.0f;
	}

	//vertex internal data
	int size; ///< Size of the vertex (number of elements). 
	int sizeFrequency; ///< Size of the vertex in the frequency subspace (number of elements). 

	float *_data; ///< Vertex data. 
	int _memoryStep[4]; ///< Memory stride width. 
};