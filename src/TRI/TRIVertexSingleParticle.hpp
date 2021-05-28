/**
 * @file TRIVertexSingleParticle.hpp
 * @author Finn Lasse Buessen
 * @brief Single-particle vertex implementation for time reversal invariant models.
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <istream>
#include "lib/Assert.hpp"
#include "FrgCommon.hpp"

/**
 * @brief Single-particle vertex implementation for time reversal invariant models.
 */
struct TRIVertexSingleParticle
{
public:
	/**
	 * @brief Construct a new TRIVertexSingleParticle object and initialize all values to zero. 
	 */
	TRIVertexSingleParticle()
	{
		//set size
		size = FrgCommon::frequency().size;

		//alloc and init memory
		_data = new float[size];
		for (int i = 0; i < size; ++i) _data[i] = 0.0f;
	}

	/**
	 * @brief Destroy the TRIVertexSingleParticle object. 
	 */
	~TRIVertexSingleParticle()
	{
		delete[] _data;
	}

	/**
	 * @brief Expand a linear iterator in the range [0,TRIVertexSingleParticle::size).
	 * 
	 * @param[in] iterator Iterator to expand. 
	 * @param[out] w Frequency argument described by the iterator. 
	 */
	void expandIterator(const int iterator, float &w) const
	{
		ASSERT(iterator >= 0 && iterator < size);

		w = FrgCommon::frequency()._data[iterator];
	}

	/**
	 * @brief Directly access a vertex value by reference via a linear iterator in the range [0,TRIVertexSingleParticle::size). 
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
	 * @brief Access vertex value at arbitrary frequency value by performing a linear interpolation on the FrequencyDiscretization. 
	 * 
	 * @param w Frequency argument. 
	 * @return float Vertex value. 
	 */
	float getValue(float w) const
	{
		int lower, upper;
		float bias;
		float sign = 1.0f;

		if (w < 0)
		{
			w = -w;
			sign = -1.0f;
		}

		FrgCommon::frequency().interpolateOffset(w, lower, upper, bias);
		return sign * ((1 - bias) * _directAccess(lower) + bias * _directAccess(upper));
	}

	/**
	 * @brief Access vertex value at given frequency mesh point. 
	 * 
	 * @param wOffset Linear offset on the frequency mesh. 
	 * @return float& Vertex value. 
	 */
	float &_directAccess(const int wOffset) const
	{
		ASSERT(wOffset >= 0 && wOffset < FrgCommon::frequency().size);

		return _data[wOffset];
	}

	int size; ///< Total number of vertex elements. 
	float *_data; ///< Vertex data. 
};