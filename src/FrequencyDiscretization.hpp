/**
 * @file FrequencyDiscretization.hpp
 * @author Finn Lasse Buessen
 * @brief Discretization of Matsubara frequency space. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "lib/Exception.hpp"
#include "lib/Log.hpp"
#include "lib/Assert.hpp"

#pragma region FrequencyIterator
/**
 * @brief Frequency iterator. 
 */
struct FrequencyIterator
{
	/**
	 * @brief Construct a new FrequencyIterator object, and initialize to a specific frequency value. 
	 * 
	 * @param p Pointer to the initial value. 
	 */
	FrequencyIterator(float *p)
	{
		_pointer = p;
	}

	/**
	 * @brief Dereference operator. 
	 * 
	 * @return float Value of the frequency which the iterator points to. 
	 */
	float operator*() const
	{
		return *_pointer;
	}

	/**
	 * @brief Iterator comparison. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Returns true if the iterators point to the same frequency value, otherwise returns false. 
	 */
	bool operator==(const FrequencyIterator &rhs) const
	{
		return _pointer == rhs._pointer;
	}

	/**
	 * @brief Negative iterator comparison. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Returns false if the iterators point to the same frequency value, otherwise returns true. 
	 */
	bool operator!=(const FrequencyIterator &rhs) const
	{
		return _pointer != rhs._pointer;
	}

	/**
	 * @brief Prefix increment operator. 
	 * 
	 * @return FrequencyIterator& Reference to self. 
	 */
	FrequencyIterator& operator++()
	{
		++_pointer;
		return *this;
	}

	/**
	 * @brief Prefix decrement operator. 
	 * 
	 * @return FrequencyIterator& Reference to self. 
	 */
	FrequencyIterator& operator--()
	{
		--_pointer;
		return *this;
	}
	
	/**
	 * @brief Iterator addition operator. 
	 * 
	 * @param rhs Number of steps to increment iterator. 
	 * @return FrequencyIterator Iterator incremented by the specified number of steps. 
	 */
	FrequencyIterator operator+(const int rhs) const
	{
		return FrequencyIterator(_pointer + rhs);
	}

	/**
	 * @brief Iterator subtraction operator. 
	 * 
	 * @param rhs Number of steps to decrement iterator. 
	 * @return FrequencyIterator Iterator decremented by the specified number of steps. 
	 */
	FrequencyIterator operator-(const int rhs) const
	{
		return FrequencyIterator(_pointer - rhs);
	}

	/**
	 * @brief Greater comparison operator. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Return true if the iterator is greater than the specified value, otherwise return false. 
	 */
	bool operator>(const FrequencyIterator &rhs) const
	{
		return _pointer > rhs._pointer;
	}
	
	/**
	 * @brief Greater or equal comparison operator. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Return true if the iterator is greater or equal than the specified value, otherwise return false. 
	 */
	bool operator>=(const FrequencyIterator &rhs) const
	{
		return _pointer >= rhs._pointer;
	}
	
	/**
	 * @brief Lesser comparison operator. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Return true if the iterator is lesser than the specified value, otherwise return false. 
	 */

	bool operator<(const FrequencyIterator &rhs) const
	{
		return _pointer < rhs._pointer;
	}
	
	/**
	 * @brief Lesser or equal comparison operator. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Return true if the iterator is lesser or equal than the specified value. 
	 */
	bool operator<=(const FrequencyIterator &rhs) const
	{
		return _pointer <= rhs._pointer;
	}

	float *_pointer; ///< Pointer which the iterator currently points to. 
};
#pragma endregion

/**
 * @brief Discretization of Matsubara frequency space. 
 * @details This structure represents a discretization of Matsubara frequency space base on a list of specified mesh points. 
 * It provides methods to iterate the frequency mesh, to find closest mesh points, and to linearly interpolate arbitrary values inbetween the mesh points. 
 * The mesh is mirror symmetric around the origin.  
 */
struct FrequencyDiscretization
{
public:
	/**
	 * @brief Construct a frequency discretization based on a list of specified mesh points. 
	 * The list of specified mesh points must be in ascending order and positive definite. 
	 * The symmetry-related negative values are automatically generated. 
	 * 
	 * @param values List of mesh points. 
	 */
	FrequencyDiscretization(const std::vector<float> &values)
	{
		//Ensure that the frequency values are in ascending order
		ASSERT(std::is_sorted(values.begin(), values.end()));

		//Ensure that discretization contains sufficiently many frequencies
		if (values.size() < 2) throw Exception(Exception::Type::ArgumentError, "FrequencyDiscretization must contain at least two frequency values");
		
		//Allocate memory and store frequencies
		size = int(values.size());
		_dataNegative = new float[2 * size];
		_data = _dataNegative + size;
		
		for (int i = 0; i < size; ++i)
		{
			ASSERT(values[i] > 0);
			_data[-i- 1] = -values[i];
			_data[i] = values[i];
		}

		Log::log << Log::LogLevel::Debug<< "Initialized frequency grid with mesh values" << Log::endl;
		for (auto i = beginNegative(); i != end(); ++i)	Log::log << "\t" << *i << Log::endl;
	}

	/**
	 * @brief Destroy the FrequencyDiscretization object. 
	 */
	~FrequencyDiscretization()
	{
		delete[] _dataNegative;
	}

	/**
	 * @brief Retrieve iterator to the first positive mesh point; This is the positive value with the smallest absolute value. 
	 * 
	 * @return FrequencyIterator Iterator to the first positive mesh point. 
	 */
	FrequencyIterator begin() const
	{
		return FrequencyIterator(_data);
	}

	/**
	 * @brief Retrieve iterator to the first negative mesh point; This is the negative value with the largest absolute value. 
	 * 
	 * @return FrequencyIterator Iterator to the first negative mesh point.
	 */
	FrequencyIterator beginNegative() const
	{
		return FrequencyIterator(_dataNegative);
	}

	/**
	 * @brief Retrieve iterator to the last mesh point; This is the positive value with the largest absolute value. 
	 * 
	 * @return FrequencyIterator Iterator to the last mesh point. 
	 */
	FrequencyIterator last() const
	{
		return FrequencyIterator(_data + size - 1);
	}

	/**
	 * @brief Retrieve iterator to the last+1 mesh point. 
	 * 
	 * @return FrequencyIterator Iterator to the last+1 mesh point. 
	 */
	FrequencyIterator end() const
	{
		return FrequencyIterator(_data + size);
	}

	/**
	 * @brief Retrieve an iterator to the closest mesh point that is lesser than the specified frequency value. 
	 * If no lesser mesh point exists, returns an iterator to the closest mesh point.  
	 * 
	 * @param w Designated upper frequency bound. 
	 * @return FrequencyIterator Iterator to the closest mesh point. 
	 */
	FrequencyIterator lesser(const float w) const
	{
		if (w < 0)
		{
			FrequencyIterator it = greater(-w);
			it._pointer = _data - (it._pointer - _data + 1);
			return it;
		}
		else
		{
			if (w <= _data[0]) return FrequencyIterator(_data);
			for (int i = 1; i < size; ++i)
			{
				if (_data[i] > w) return FrequencyIterator(_data + i - 1);
			}
			return FrequencyIterator(_data + size - 1);
		}
	}

	/**
	 * @brief Retrieve an iterator to the closest mesh point that is greater than the specified frequency value. 
	 * If no greater mesh point exists, returns an iterator to the closest mesh point.  
	 * 
	 * @param w Designated upper frequency bound. 
	 * @return FrequencyIterator Iterator to the closest mesh point. 
	 */
	FrequencyIterator greater(const float w) const
	{
		if (w < 0)
		{
			FrequencyIterator it = lesser(-w);
			it._pointer = _data - (it._pointer - _data + 1);
			return it;
		}
		else
		{
			if (w <= _data[0]) return FrequencyIterator(_data);
			for (int i = 1; i < size; ++i)
			{
				if (_data[i] > w) return FrequencyIterator(_data + i);
			}
			return FrequencyIterator(_data + size - 1);
		}
	}

	/**
	 * @brief Return the number of iterator increments of a mesh point associated with a given frequency value, relative to the first positive mesh point. 
	 * The value of the specified frequency must be positive. 
	 * If no mesh point with that value exists, returns the index of the last mesh point. 
	 * 
	 * @param w Frequency value. Must be positive. 
	 * @return int Index of the mesh point. 
	 */
	int offset(const float w) const
	{
		ASSERT(w >= 0);

		if (w <= _data[0]) return 0;
		for (int i = 1; i < size; ++i)
		{
			if (_data[i] >= w) return i;
		}
		return size - 1;
	}	

	/**
	 * @brief Perform an interpolation between mesh points for an arbitrary positive frequency. 
	 * 
	 * @param[in] w Frequency value to interpolate to. Must be positive. 
	 * @param[out] lowerOffset Number of iterator increments of the lesser mesh point, relative to the first positive mesh point. 
	 * @param[out] upperOffset Number of iterator increments of the greater mesh point, relative to the first positive mesh point. 
	 * @param[out] bias Linear interpolation weight. Zero, if `w` matches the mesh point at `lowerOffset`. One. if `w` matches the mesh point at `upperOffset`. 
	 */
	void interpolateOffset(const float w, int &lowerOffset, int &upperOffset, float &bias) const
	{
		ASSERT(w >= 0);
		ASSERT(&lowerOffset != &upperOffset)

		if (w <= _data[0])
		{
			lowerOffset = 0; 
			upperOffset = 0;
			bias = 0.0f;
			return;
		}
		for (int i = 1; i < size; ++i)
		{
			if (_data[i] > w)
			{
				upperOffset = i;
				lowerOffset = i - 1;
				bias = (w - _data[lowerOffset]) / (_data[upperOffset] - _data[lowerOffset]);
				return;
			}
		}
		lowerOffset = size - 1;
		upperOffset = size - 1;
		bias = 0.0f;
		return;

		ASSERT(lowerOffset >= 0 && lowerOffset < size);
		ASSERT(upperOffset >= 0 && upperOffset < size);
		ASSERT(bias >= 0.0f && bias <= 1.0f);
	}

	int size; ///< Number of positive mesh points. 
	float *_data; ///< Pointer to the first positive mesh point. Stored contiuously after FrequencyDiscretization::_dataNegative. 
	float *_dataNegative; ///< Pointer to the first negative mesh point. 
};
