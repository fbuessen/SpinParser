/**
 * @file SU2VertexTwoParticle.hpp
 * @author Finn Lasse Buessen
 * @brief Two-particle vertex implementation for SU(2) models.
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
 * Each support value is also complemented by a sign factor that is either +1 or -1 and by the information whether a lattice site exchange should be performed upon accessing the support value. 
 * 
 * @tparam size Number of support sites for the interpolation. 
 */
template <int size> struct SU2VertexTwoParticleAccessBuffer
{
public:
	/**
	 * @brief Construct an uninitialized SU2VertexTwoParticleAccessBuffer object. 
	 */
	SU2VertexTwoParticleAccessBuffer() : siteExchange(false) {}

	int frequencyOffsets[size]; ///< Linear memory offset (number of elements) in the frequency dimensions of the two-particle vertex.  
	float frequencyWeights[size]; ///< Weight factors of the support values. 
	int signFlag[size]; ///< Sign factors of the support values. 
	bool siteExchange; ///< Site exchange indicator. 
};

/**
 * @brief Two-particle vertex implementation for SU(2) models.
 */
struct SU2VertexTwoParticle
{
public:
	/**
	 * @brief Enumeration of the different vertex channels. 
	 */
	enum struct Symmetry : int
	{
		Spin = 0, ///< Spin-type vertex. 
		Density = 1 ///< Density-type vertex. 
	};

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
	 * @brief Construct a new SU2VertexTwoParticle object and initialize all entries to zero. 
	 */
	SU2VertexTwoParticle()
	{
		//store width in all memory dimensions
		_memoryStepLattice = FrgCommon::lattice().size;
		_memoryStepLatticeT = _memoryStepLattice * FrgCommon::frequency().size;

		sizeFrequency = FrgCommon::frequency().size * FrgCommon::frequency().size * (FrgCommon::frequency().size + 1) / 2;
		size = FrgCommon::lattice().size * sizeFrequency;

		//alloc and init memory
		_dataSS = new float[size];
		_dataDD = new float[size];
		memset(_dataSS, 0, sizeof(float) * size);
		memset(_dataDD, 0, sizeof(float) * size);
	}

	/**
	 * @brief Destroy the SU2VertexTwoParticle object. 
	 */
	~SU2VertexTwoParticle()
	{
		delete[] _dataSS;
		delete[] _dataDD;
	}

	/**
	 * @brief Expand a linear iterator in the range [0,size) that iterates over all vertex entries. 
	 * 
	 * @param[in] iterator Linear iterator. 
	 * @param[out] i1 Lattice site iterator. 
	 * @param[out] s First frequency argument. 
	 * @param[out] t Second frequency argument. 
	 * @param[out] u Third frequency argument. 
	 */
	void expandIterator(int iterator, LatticeIterator &i1, float &s, float &t, float &u) const
	{
		ASSERT(iterator >= 0 && iterator < size);
		ASSERT(&s != &t);
		ASSERT(&t != &u);
		ASSERT(&s != &u);

		int su = iterator / _memoryStepLatticeT;
		iterator = iterator % _memoryStepLatticeT;
		t = FrgCommon::frequency()._data[iterator / _memoryStepLattice];
		i1 = FrgCommon::lattice().fromParametrization(iterator % _memoryStepLattice);

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
	 * @param symmetry Vertex channel. 
	 * @return float& Vertex value. 
	 */
	float &getValueRef(const int iterator, const SU2VertexTwoParticle::Symmetry symmetry) const
	{
		if (symmetry == SU2VertexTwoParticle::Symmetry::Spin) return _dataSS[iterator];
		else return _dataDD[iterator];
	}

	/**
	 * @brief Access vertex value at arbitrary lattice sites, frequencies, and symmetry. 
	 * 
	 * @param i1 First lattice site argument. 
	 * @param i2 Second lattice site argument. 
	 * @param s First frequency argument. 
	 * @param t Second frequency argument. 
	 * @param u Third frequency argument. 
	 * @param symmetry Vertex channel. 
	 * @param channel Frequency channel. 
	 * @return float Vertex value. 
	 */
	float getValue(LatticeIterator i1, LatticeIterator i2, float s, float t, float u, const SU2VertexTwoParticle::Symmetry symmetry, const SU2VertexTwoParticle::FrequencyChannel channel) const
	{
		//map to positive frequency sector
		if (s < 0 && u < 0)
		{
			s = -s;
			u = -u;
		}
		else
		{
			if (s < 0)
			{
				s = -s;
				LatticeIterator i3 = i1;
				i1 = i2;
				i2 = i3;
			}
			else if (u < 0)
			{
				u = -u;
				LatticeIterator i3 = i1;
				i1 = i2;
				i2 = i3;
			}
		}
		if (t < 0)
		{
			t = -t;
		}

		int siteOffset = FrgCommon::lattice().symmetryTransform(i1, i2);
		if (channel == SU2VertexTwoParticle::FrequencyChannel::S)
		{
			int exactS = FrgCommon::frequency().offset(s);

			int lowerT, upperT;
			float biasT;
			int lowerU, upperU;
			float biasU;

			FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);
			FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

			return (1 - biasU) * (
				(1 - biasT) * (_directAccessMapFrequencyExchange(siteOffset, exactS, lowerT, lowerU, symmetry)) + biasT * (_directAccessMapFrequencyExchange(siteOffset, exactS, upperT, lowerU, symmetry))
				) + biasU * (
				(1 - biasT) * (_directAccessMapFrequencyExchange(siteOffset, exactS, lowerT, upperU, symmetry)) + biasT * (_directAccessMapFrequencyExchange(siteOffset, exactS, upperT, upperU, symmetry))
					);
		}
		else if (channel == SU2VertexTwoParticle::FrequencyChannel::T)
		{
			int exactT = FrgCommon::frequency().offset(t);

			int lowerS, upperS;
			float biasS;
			int lowerU, upperU;
			float biasU;

			FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
			FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

			return (1 - biasU) * (
				(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, exactT, lowerU, symmetry)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, exactT, lowerU, symmetry))
				) + biasU * (
				(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, exactT, upperU, symmetry)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, exactT, upperU, symmetry))
					);
		}
		else if (channel == SU2VertexTwoParticle::FrequencyChannel::U)
		{
			int exactU = FrgCommon::frequency().offset(u);

			int lowerS, upperS;
			float biasS;
			int lowerT, upperT;
			float biasT;

			FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
			FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);

			return (1 - biasT) * (
				(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, lowerT, exactU, symmetry)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, lowerT, exactU, symmetry))
				) + biasT * (
				(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, upperT, exactU, symmetry)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, upperT, exactU, symmetry))
					);
		}
		else if (channel == SU2VertexTwoParticle::FrequencyChannel::None)
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

			return
				(1 - biasU) * (
				(1 - biasT) * (
					(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, lowerT, lowerU, symmetry)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, lowerT, lowerU, symmetry))
					) + biasT * (
					(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, upperT, lowerU, symmetry)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, upperT, lowerU, symmetry))
						)
					) + biasU * (
					(1 - biasT) * (
						(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, lowerT, upperU, symmetry)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, lowerT, upperU, symmetry))
						) + biasT * (
						(1 - biasS) * (_directAccessMapFrequencyExchange(siteOffset, lowerS, upperT, upperU, symmetry)) + biasS * (_directAccessMapFrequencyExchange(siteOffset, upperS, upperT, upperU, symmetry))
							)
						);
		}
		else if (channel == SU2VertexTwoParticle::FrequencyChannel::All)
		{
			int exactS = FrgCommon::frequency().offset(s);
			int exactT = FrgCommon::frequency().offset(t);
			int exactU = FrgCommon::frequency().offset(u);
			return _directAccessMapFrequencyExchange(siteOffset, exactS, exactT, exactU, symmetry);
		}
		else
		{
			throw Exception(Exception::Type::ArgumentError, "Invalid interpolation channel");
			return 0.0f;
		}
	}

	/**
	 * @brief Access vertex value at arbitrary lattice sites and symmetry via a given access buffer. 
	 * 
	 * @tparam n Number of support sites in the access buffer. 
	 * @param i1 First lattice site argument. 
	 * @param i2 Second lattice site argument. 
	 * @param symmetry Vertex channel. 
	 * @param accessBuffer Access buffer. 
	 * @return float Vertex value. 
	 */
	template <int n> float getValue(const LatticeIterator i1, const LatticeIterator i2, const SU2VertexTwoParticle::Symmetry symmetry, const SU2VertexTwoParticleAccessBuffer<n> &accessBuffer) const
	{
		int siteOffset = (accessBuffer.siteExchange) ? FrgCommon::lattice().symmetryTransform(i2, i1) : FrgCommon::lattice().symmetryTransform(i1, i2);

		float value = 0.0f;
		if (symmetry == SU2VertexTwoParticle::Symmetry::Spin)
		{
			for (int i = 0; i < n; ++i) value += accessBuffer.frequencyWeights[i] * _dataSS[accessBuffer.frequencyOffsets[i] + siteOffset];
		}
		else
		{
			for (int i = 0; i < n; ++i) value += accessBuffer.signFlag[i] * accessBuffer.frequencyWeights[i] * _dataDD[accessBuffer.frequencyOffsets[i] + siteOffset];
		}

		return value;
	}

	/**
	 * @brief Access vertex value locally at arbitrary symmetry via a given access buffer
	 * 
	 * @tparam n Number of support sites in the access buffer. 
	 * @param symmetry Vertex channel. 
	 * @param accessBuffer Access buffer. 
	 * @return float Vertex value. 
	 */
	template <int n> float getValueLocal(const SU2VertexTwoParticle::Symmetry symmetry, const SU2VertexTwoParticleAccessBuffer<n> &accessBuffer) const
	{
		float value = 0.0f;
		if (symmetry == SU2VertexTwoParticle::Symmetry::Spin)
		{
			for (int i = 0; i < n; ++i) value += accessBuffer.frequencyWeights[i] * _dataSS[accessBuffer.frequencyOffsets[i]];
		}
		else
		{
			for (int i = 0; i < n; ++i) value += accessBuffer.signFlag[i] * accessBuffer.frequencyWeights[i] * _dataDD[accessBuffer.frequencyOffsets[i]];
		}

		return value;
	}

	/**
	 * @brief Bundled vertex access on all lattice sites and symmetries simultaneously via a given access buffer. 
	 * 
	 * @tparam n Number of support sites in the access buffer. 
	 * @param[in] accessBuffer Access buffer. 
	 * @param[out] superbundle Vertex value bundle. 
	 */
	template <int n> void getValueSuperbundle(const SU2VertexTwoParticleAccessBuffer<n> &accessBuffer, ValueSuperbundle<float, 2> &superbundle) const
	{
		superbundle.reset();
		const LatticeSiteDescriptor *sites = (accessBuffer.siteExchange) ? FrgCommon::lattice().getInvertedSites() : FrgCommon::lattice().getSites();

		for (int i = 0; i < n; ++i)
		{
			float weight = accessBuffer.frequencyWeights[i];
			float signedWeight = accessBuffer.signFlag[i] * accessBuffer.frequencyWeights[i];
			int frequencyOffset = accessBuffer.frequencyOffsets[i];
			int size = FrgCommon::lattice().size;

			for (int j = 0; j < size; ++j)
			{
				superbundle.bundle(0)[j] += weight * _dataSS[frequencyOffset + sites[j].rid];
				superbundle.bundle(1)[j] += signedWeight * _dataDD[frequencyOffset + sites[j].rid];
			}
		}
	}

	/**
	 * @brief Generate an access buffer for a set of frequencies where one of them (specified by channel) exactly lies on the frequency mesh. 
	 * Frequency channel must be either FrequencyChannel::S, FrequencyChannel::T, FrequencyChannel::U. 
	 * 
	 * @param s First frequency argument. 
	 * @param t Second frequency argument. 
	 * @param u Third frequency argument. 
	 * @param channel Frequency channel. 
	 * @return SU2VertexTwoParticleAccessBuffer<4> Access buffer. 
	 */
	SU2VertexTwoParticleAccessBuffer<4> generateAccessBuffer(float s, float t, float u, const SU2VertexTwoParticle::FrequencyChannel channel) const
	{
		ASSERT(channel == FrequencyChannel::S || channel == FrequencyChannel::T || channel == FrequencyChannel::U);

		SU2VertexTwoParticleAccessBuffer<4> accessBuffer;

		//map to positive frequency sector
		if (s < 0 && u < 0)
		{
			s = -s;
			u = -u;
		}
		else
		{
			if (s < 0)
			{
				s = -s;
				accessBuffer.siteExchange = true;
			}
			else if (u < 0)
			{
				u = -u;
				accessBuffer.siteExchange = true;
			}
		}
		if (t < 0)
		{
			t = -t;
		}

		//interpolate frequency
		if (channel == SU2VertexTwoParticle::FrequencyChannel::S)
		{
			int exactS = FrgCommon::frequency().offset(s);
			int lowerT, upperT;
			float biasT;
			int lowerU, upperU;
			float biasU;
			FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);
			FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

			accessBuffer.frequencyWeights[0] = (1 - biasU) * (1 - biasT);
			accessBuffer.frequencyOffsets[0] = _generateAccessBufferOffset(exactS, lowerT, lowerU, accessBuffer.signFlag[0]);
			accessBuffer.frequencyWeights[1] = (1 - biasU) * biasT;
			accessBuffer.frequencyOffsets[1] = _generateAccessBufferOffset(exactS, upperT, lowerU, accessBuffer.signFlag[1]);
			accessBuffer.frequencyWeights[2] = biasU * (1 - biasT);
			accessBuffer.frequencyOffsets[2] = _generateAccessBufferOffset(exactS, lowerT, upperU, accessBuffer.signFlag[2]);
			accessBuffer.frequencyWeights[3] = biasU * biasT;
			accessBuffer.frequencyOffsets[3] = _generateAccessBufferOffset(exactS, upperT, upperU, accessBuffer.signFlag[3]);
		}
		else if (channel == SU2VertexTwoParticle::FrequencyChannel::T)
		{
			int exactT = FrgCommon::frequency().offset(t);
			int lowerS, upperS;
			float biasS;
			int lowerU, upperU;
			float biasU;
			FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
			FrgCommon::frequency().interpolateOffset(u, lowerU, upperU, biasU);

			accessBuffer.frequencyWeights[0] = (1 - biasU) * (1 - biasS);
			accessBuffer.frequencyOffsets[0] = _generateAccessBufferOffset(lowerS, exactT, lowerU, accessBuffer.signFlag[0]);
			accessBuffer.frequencyWeights[1] = (1 - biasU) * biasS;
			accessBuffer.frequencyOffsets[1] = _generateAccessBufferOffset(upperS, exactT, lowerU, accessBuffer.signFlag[1]);
			accessBuffer.frequencyWeights[2] = biasU * (1 - biasS);
			accessBuffer.frequencyOffsets[2] = _generateAccessBufferOffset(lowerS, exactT, upperU, accessBuffer.signFlag[2]);
			accessBuffer.frequencyWeights[3] = biasU * biasS;
			accessBuffer.frequencyOffsets[3] = _generateAccessBufferOffset(upperS, exactT, upperU, accessBuffer.signFlag[3]);

		}
		else if (channel == SU2VertexTwoParticle::FrequencyChannel::U)
		{
			int exactU = FrgCommon::frequency().offset(u);
			int lowerS, upperS;
			float biasS;
			int lowerT, upperT;
			float biasT;
			FrgCommon::frequency().interpolateOffset(s, lowerS, upperS, biasS);
			FrgCommon::frequency().interpolateOffset(t, lowerT, upperT, biasT);

			accessBuffer.frequencyWeights[0] = (1 - biasT) * (1 - biasS);
			accessBuffer.frequencyOffsets[0] = _generateAccessBufferOffset(lowerS, lowerT, exactU, accessBuffer.signFlag[0]);
			accessBuffer.frequencyWeights[1] = (1 - biasT) * biasS;
			accessBuffer.frequencyOffsets[1] = _generateAccessBufferOffset(upperS, lowerT, exactU, accessBuffer.signFlag[1]);
			accessBuffer.frequencyWeights[2] = biasT * (1 - biasS);
			accessBuffer.frequencyOffsets[2] = _generateAccessBufferOffset(lowerS, upperT, exactU, accessBuffer.signFlag[2]);
			accessBuffer.frequencyWeights[3] = biasT * biasS;
			accessBuffer.frequencyOffsets[3] = _generateAccessBufferOffset(upperS, upperT, exactU, accessBuffer.signFlag[3]);
		}

		return accessBuffer;
	}

	/**
	 * @brief Generate an access buffer for an arbitrary set of frequencies.  
	 * 
	 * @param s First frequency argument. 
	 * @param t Second frequency argument. 
	 * @param u Third frequency argument. 
	 * @return SU2VertexTwoParticleAccessBuffer<8> Access buffer. 
	 */
	SU2VertexTwoParticleAccessBuffer<8> generateAccessBuffer(float s, float t, float u) const
	{
		SU2VertexTwoParticleAccessBuffer<8> accessBuffer;

		//map to positive frequency sector
		if (s < 0 && u < 0)
		{
			s = -s;
			u = -u;
		}
		else
		{
			if (s < 0)
			{
				s = -s;
				accessBuffer.siteExchange = true;
			}
			else if (u < 0)
			{
				u = -u;
				accessBuffer.siteExchange = true;
			}
		}
		if (t < 0)
		{
			t = -t;
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

		accessBuffer.frequencyWeights[0] = (1 - biasT) * (1 - biasS) * (1 - biasU);
		accessBuffer.frequencyOffsets[0] = _generateAccessBufferOffset(lowerS, lowerT, lowerU, accessBuffer.signFlag[0]);
		accessBuffer.frequencyWeights[1] = (1 - biasT) * biasS * (1 - biasU);
		accessBuffer.frequencyOffsets[1] = _generateAccessBufferOffset(upperS, lowerT, lowerU, accessBuffer.signFlag[1]);
		accessBuffer.frequencyWeights[2] = biasT * (1 - biasS) * (1 - biasU);
		accessBuffer.frequencyOffsets[2] = _generateAccessBufferOffset(lowerS, upperT, lowerU, accessBuffer.signFlag[2]);
		accessBuffer.frequencyWeights[3] = biasT * biasS * (1 - biasU);
		accessBuffer.frequencyOffsets[3] = _generateAccessBufferOffset(upperS, upperT, lowerU, accessBuffer.signFlag[3]);
		accessBuffer.frequencyWeights[4] = (1 - biasT) * (1 - biasS) * biasU;
		accessBuffer.frequencyOffsets[4] = _generateAccessBufferOffset(lowerS, lowerT, upperU, accessBuffer.signFlag[4]);
		accessBuffer.frequencyWeights[5] = (1 - biasT) * biasS * biasU;
		accessBuffer.frequencyOffsets[5] = _generateAccessBufferOffset(upperS, lowerT, upperU, accessBuffer.signFlag[5]);
		accessBuffer.frequencyWeights[6] = biasT * (1 - biasS) * biasU;
		accessBuffer.frequencyOffsets[6] = _generateAccessBufferOffset(lowerS, upperT, upperU, accessBuffer.signFlag[6]);
		accessBuffer.frequencyWeights[7] = biasT * biasS * biasU;
		accessBuffer.frequencyOffsets[7] = _generateAccessBufferOffset(upperS, upperT, upperU, accessBuffer.signFlag[7]);

		return accessBuffer;
	}

	/**
	 * @brief Directly access a vertex via given frequency and site offsets, where sOffset may be lesser than uOffset. 
	 * 
	 * @param siteOffset Lattice site offset (number of elements). 
	 * @param sOffset First frequency offset (number of elements). 
	 * @param tOffset Second frequency offset (number of elements). 
	 * @param uOffset Third frequency offset (number of elements). 
	 * @param symmetry Vertex channel. 
	 * @return float Vertex value. 
	 */
	float _directAccessMapFrequencyExchange(const int siteOffset, const int sOffset, const int tOffset, const int uOffset, const SU2VertexTwoParticle::Symmetry symmetry) const
	{
		ASSERT(siteOffset >= 0 && siteOffset < FrgCommon::lattice().size);
		ASSERT(sOffset >= 0 && sOffset < FrgCommon::frequency().size);
		ASSERT(tOffset >= 0 && tOffset < FrgCommon::frequency().size);
		ASSERT(uOffset >= 0 && uOffset < FrgCommon::frequency().size);


		if (sOffset >= uOffset) return _directAccess(siteOffset, sOffset, tOffset, uOffset, symmetry);
		else
		{
			if (symmetry == SU2VertexTwoParticle::Symmetry::Spin) return _directAccess(siteOffset, uOffset, tOffset, sOffset, symmetry);
			else return -_directAccess(siteOffset, uOffset, tOffset, sOffset, symmetry);
		}
	}

	/**
	 * @brief Directly access a vertex via given frequency and site offsets, where sOffset >= uOffset. 
	 * 
	 * @param siteOffset Lattice site offset (number of elements). 
	 * @param sOffset First frequency offset (number of elements). 
	 * @param tOffset Second frequency offset (number of elements). 
	 * @param uOffset Third frequency offset (number of elements). 
	 * @param symmetry Vertex channel. 
	 * @return float Vertex value. 
	 */
	float _directAccess(const int siteOffset, const int sOffset, const int tOffset, const int uOffset, const SU2VertexTwoParticle::Symmetry symmetry) const
	{
		ASSERT(siteOffset >= 0 && siteOffset < FrgCommon::lattice().size);
		ASSERT(sOffset >= 0 && sOffset < FrgCommon::frequency().size);
		ASSERT(tOffset >= 0 && tOffset < FrgCommon::frequency().size);
		ASSERT(uOffset >= 0 && uOffset < FrgCommon::frequency().size);
		ASSERT(sOffset >= uOffset);

		if (symmetry == SU2VertexTwoParticle::Symmetry::Spin) return _dataSS[_memoryStepLatticeT * (sOffset * (sOffset + 1) / 2 + uOffset) + _memoryStepLattice * tOffset + siteOffset];
		else return _dataDD[_memoryStepLatticeT * (sOffset * (sOffset + 1) / 2 + uOffset) + _memoryStepLattice * tOffset + siteOffset];
	}

	/**
	 * @brief Calculate the total memory offset (number of elements) from given frequency offsets, where sOffset may be lesser than uOffset. 
	 * 
	 * @param[in] sOffset First frequency offset (number of elements). 
	 * @param[in] tOffset Second frequency offset (number of elements). 
	 * @param[in] uOffset Third frequency offset (number of elements). 
	 * @param[out] signFlag Sign flag to be stored in the access buffer. 
	 * @return int Memory offset (number of elements). 
	 */
	int _generateAccessBufferOffset(const int sOffset, const int tOffset, const int uOffset, int &signFlag) const
	{
		ASSERT(sOffset >= 0 && sOffset < FrgCommon::frequency().size);
		ASSERT(tOffset >= 0 && tOffset < FrgCommon::frequency().size);
		ASSERT(uOffset >= 0 && uOffset < FrgCommon::frequency().size);

		if (sOffset < uOffset)
		{
			signFlag = -1;
			return _memoryStepLatticeT * (uOffset * (uOffset + 1) / 2 + sOffset) + tOffset * _memoryStepLattice;
		}
		else
		{
			signFlag = 1;
			return _memoryStepLatticeT * (sOffset * (sOffset + 1) / 2 + uOffset) + tOffset * _memoryStepLattice;
		}
	}

	int size; ///< Size of the vertex per vertex channel (number of elements). 
	int sizeFrequency; ///< Size of the vertex per vertex channel in the frequency subspace (number of elements). 

	float *_dataSS; ///< Spin channel of the vertex. 
	float *_dataDD; ///< Density channel of the vertex. 
	int _memoryStepLatticeT; ///< Memory stride width in the last-2 dimension. 
	int _memoryStepLattice; ///< Memory stride width in the last-1 dimension. 
};