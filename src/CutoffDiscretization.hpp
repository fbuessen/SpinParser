/**
 * @file CutoffDiscretization.hpp
 * @author Finn Lasse Buessen
 * @brief Representation of a discretized frequency cutoff axis. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <vector>
#include "lib/Exception.hpp"

#pragma region CutoffIterator
/**
 * @brief Iterator over discrete cutoff values. 
 */
struct CutoffIterator
{
public:
	/**
	 * @brief Construct a new CutoffIterator object pointing to specific value. 
	 * 
	 * @param p Pointer to initial value. 
	 */
	CutoffIterator(float* p) : _pointer(p) {};

	/**
	 * @brief Iterator dereference. 
	 * 
	 * @return float Cutoff value. 
	 */
	float operator*() const
	{
		return *_pointer;
	}

	/**
	 * @brief Iterator comparison
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Returns true of iterators are equal, returns false otherwise.  
	 */
	bool operator==(const CutoffIterator& rhs) const
	{
		return _pointer == rhs._pointer;
	}
	
	/**
	 * @brief Negative iterator comparison
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Returns true of iterators are unequal, returns false otherwise.   
	 */
	bool operator!=(const CutoffIterator& rhs) const
	{
		return _pointer != rhs._pointer;
	}
	
	/**
	 * @brief Prefix increment. 
	 * 
	 * @return CutoffIterator& Reference to self. 
	 */
	CutoffIterator& operator++()
	{
		++_pointer;
		return *this;
	}

private:
	float *_pointer; ///< Pointer which the iterator currently points to. 
};
#pragma endregion

/**
 * @brief Representation of frequency axis cutoff discretization. 
 */
struct CutoffDiscretization
{
public:
	/**
	 * @brief Construct a new CutoffDiscretization object from a list of cutoff values. 
	 * 
	 * @param values List of cutoff values to use for discretization. 
	 */
	CutoffDiscretization(const std::vector<float> &values)
	{
		//Ensure that discretization contains sufficiently many cutoff values
		if (values.size() < 2) throw Exception(Exception::Type::ArgumentError, "CutoffDiscretization must contain at least two frequency values");

		_size = int(values.size());
		_data = new float[values.size()];
		memcpy(_data, values.data(), values.size() * sizeof(float));
	}

	/**
	 * @brief Destroy the CutoffDiscretization object
 	*/
	~CutoffDiscretization()
	{
		delete[] _data;
	}

	/**
	 * @brief Retrieve iterator to first discretization value. 
	 * 
	 * @return CutoffIterator Iterator to first discretization value. 
	 */
	CutoffIterator begin() const
	{
		return CutoffIterator(_data);
	}

	/**
	 * @brief Retrieve iterator to last discretization value. 
	 * 
	 * @return CutoffIterator Iterator to last discretization value. 
	 */
	CutoffIterator last() const
	{
		return CutoffIterator(_data + _size - 1);
	}

	/**
	 * @brief Retrieve iterator to last+1 discretization value. 
	 * 
	 * @return CutoffIterator Iterator to last+1 discretization value. 
	 */
	CutoffIterator end() const
	{
		return CutoffIterator(_data + _size);
	}

	/**
	 * @brief Retrieve iterator pointing to a specific cutoff value. 
	 * 
	 * @param cutoff Cutoff search value. 
	 * @return CutoffIterator Iterator pointing to specified cutoff value it it exists; otherwise points to last+1. 
	 */
	CutoffIterator find(const float cutoff) const
	{
		for (int i = 0; i < _size; ++i)
		{
			if (_data[i] == cutoff) return CutoffIterator(_data + i);
		}
		return end();
	}

private:
	int _size; ///< Number of cutoff values in the discretization. 
	float *_data; ///< Internal storage for discretization values. 
};