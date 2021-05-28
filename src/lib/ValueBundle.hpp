/**
 * @file ValueBundle.hpp
 * @author Finn Lasse Buessen
 * @brief Lightweight library for value arrays and collections thereof. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <cstring>

/**
 * @brief Value array implementation. The object does not hold ownership of its memory. 
 * 
 * @tparam T Underlying fundamental data type. 
 */
template <class T> struct ValueBundle
{
public:
	/**
	 * @brief Construct an empty ValueBundle object. 
	 */
	ValueBundle() : _data(nullptr), _size(0) {}

	/**
	 * @brief Construct a new ValueBundle object.
	 * 
	 * @param data Storage memory, size should at least be size * sizeof(T). Memory is not deleted on destruction of the ValueBundle.
	 * @param size Number of elemenets in the ValueBundle. 
	 */
	ValueBundle(T *data, const int size) : _data(data), _size(size) {}
	
	/**
	 * @brief Assignment operator. 
	 * 
	 * @param rhs Right hand size operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &operator=(const ValueBundle<T> &rhs)
	{
		_size = rhs._size;
		_data = rhs._data;
		return *this;
	}

	/**
	 * @brief Access the nth value stored in the value bundle. 
	 * 
	 * @param n Element number. 
	 * @return T& Reference to nth element.
	 */
	T &operator[](const int n)
	{
		return _data[n];
	}

	/**
	 * @brief Retrieve the number of elements in the value bundle. 
	 * 
	 * @return int Number of elements in the value bundle.
	 */
	int size() const
	{
		return _size;
	}

	/**
	 * @brief Retrieve the pointer to the first value stored in the bundle. 
	 * 
	 * @return const T* Pointer to the first vlaue stored in the bundle.
	 */
	T *data() const
	{
		return _data;
	}

	/**
	 * @brief Fused multiply-add, equivalent to operator+=(rhs1 * rhs2). 
	 * 
	 * @param rhs1 Scalar operand. 
	 * @param rhs2 Array operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &multAdd(const T &rhs1, const ValueBundle<T> &rhs2)
	{
		for (int i = 0; i < _size; ++i) _data[i] += rhs1 * rhs2._data[i];
		return *this;
	}

	/**
	 * @brief Fused multiply-add for scalar right-multiplication, equivalent to operator+=(rhs1 * rhs2). 
	 * 
	 * @param rhs1 Array operand. 
	 * @param rhs2 Scalar operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &multAdd(const ValueBundle<T> &rhs1, const T &rhs2)
	{
		return multAdd(rhs2, rhs1);
	}

	/**
	 * @brief Elementwise fused multiply-add. Equivalent to elementwise operator+=(rhs1 * rhs2). 
	 * 
	 * @param rhs1 Left array operand. 
	 * @param rhs2 Right array operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &multAdd(const ValueBundle<T> &rhs1, const ValueBundle<T> &rhs2)
	{
		for (int i = 0; i < _size; ++i) _data[i] += rhs1._data[i] * rhs2._data[i];
		return *this;
	}

	/**
	 * @brief Elementwise double-multiply add. Equivalent to elementwise operator+=(rhs1 * rhs2 * rhs3). 
	 * 
	 * @param rhs1 Scalar operand. 
	 * @param rhs2 Left array operand. 
	 * @param rhs3 Right array operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &multAdd(const T &rhs1, const ValueBundle<T> &rhs2, const ValueBundle<T> &rhs3)
	{
		for (int i = 0; i < _size; ++i) _data[i] += rhs1 * rhs2._data[i] * rhs3._data[i];
		return *this;
	}

	/**
	 * @brief Fused multiply-sub, equivalent to operator-=(rhs1 * rhs2). 
	 * 
	 * @param rhs1 Scalar operand. 
	 * @param rhs2 Array operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &multSub(const T &rhs1, const ValueBundle<T> &rhs2)
	{
		for (int i = 0; i < _size; ++i) _data[i] -= rhs1 * rhs2._data[i];
		return *this;
	}

	/**
	 * @brief Fused multiply-sub for scalar right-multiplication, equivalent to operator-=(rhs1 * rhs2). 
	 * 
	 * @param rhs1 Array operand. 
	 * @param rhs2 Scalar operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &multSub(const ValueBundle<T> &rhs1, const T &rhs2)
	{
		return multSub(rhs2, rhs1);
	}

	/**
	 * @brief Elementwise fused multiply-sub. Equivalent to elementwise operator-=(rhs1 * rhs2). 
	 * 
	 * @param rhs1 Left array operand. 
	 * @param rhs2 Right array operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &multSub(const ValueBundle<T> &rhs1, const ValueBundle<T> &rhs2)
	{
		for (int i = 0; i < _size; ++i) _data[i] -= rhs1._data[i] * rhs2._data[i];
		return *this;
	}

	/**
	 * @brief Elementwise double-multiply sub. Equivalent to elementwise operator-=(rhs1 * rhs2 * rhs3). 
	 * 
	 * @param rhs1 Scalar operand. 
	 * @param rhs2 Left array operand. 
	 * @param rhs3 Right array operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &multSub(const T &rhs1, const ValueBundle<T> &rhs2, const ValueBundle<T> &rhs3)
	{
		for (int i = 0; i < _size; ++i) _data[i] -= rhs1 * rhs2._data[i] * rhs3._data[i];
		return *this;
	}

	/**
	 * @brief Addition assignment. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &operator+=(const ValueBundle &rhs)
	{
		for (int i = 0; i < _size; ++i) _data[i] += rhs._data[i];
		return *this;
	}

	/**
	 * @brief Subtraction assignment. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &operator-=(const ValueBundle &rhs)
	{
		for (int i = 0; i < _size; ++i) _data[i] -= rhs._data[i];
		return *this;
	}

	/**
	 * @brief Multiplication assignment. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &operator*=(const T &rhs)
	{
		for (int i = 0; i < _size; ++i) _data[i] *= rhs;
		return *this;
	}

	/**
	 * @brief Division assignment
	 * 
	 * @param rhs Right hand side operand. 
	 * @return ValueBundle& Reference to self. 
	 */
	ValueBundle &operator/=(const T &rhs)
	{
		for (int i = 0; i < _size; ++i) _data[i] /= rhs;
		return *this;
	}

private:
	T *_data; ///< Internal data storage. Memory is not owned by the ValueBundle. 
	int _size; ///< Number of elements in the ValueBundle. 
};

/**
 * @brief Collection of ValueBundles. 
 * 
 * @tparam T Fundamental data type of the ValueBundles. 
 * @tparam n Number of ValueBundles in the collection. 
 */
template <class T, int n> struct ValueSuperbundle
{
public:
	/**
	 * @brief Construct a new ValueSuperbundle object and allocate ValueBundles. 
	 * 
	 * @param bundleSize Number of elements in each ValueBundle. 
	 */
	ValueSuperbundle(const int bundleSize) : hasOwnership(true)
	{
		for (int i = 0; i < n; ++i) bundles[i] = ValueBundle<T>(new T[bundleSize], bundleSize);
		reset();
	}

	/**
	 * @brief Copy constructor. The copy will not have ownership of the ValueBundle memory. 
	 * 
	 * @param rhs Right hand side operand. 
	 */
	ValueSuperbundle(const ValueSuperbundle &rhs) : hasOwnership(false)
	{
		for (int i = 0; i < n; ++i) bundles[i] = ValueBundle<T>(rhs.bundles[i].data(), rhs.bundles[i].size());
	}

	/**
	 * @brief Destroy the ValueSuperbundle object
	 */
	~ValueSuperbundle()
	{
		if (hasOwnership)
		{
			for (int i = 0; i < n; ++i) delete[] bundles[i].data();
		}
	}

	/**
	 * @brief Return reference to ValueBundle. 
	 * 
	 * @param m Id of the ValueBundle. 
	 * @return ValueBundle<T>& Reference to the m-th ValueBundle. 
	 */
	ValueBundle<T> &bundle(const int m)
	{
		return bundles[m];
	}

	/**
	 * @brief Write zeros to all ValueBundles. 
	 * 
	 * @return ValueSuperbundle& Reference to self. 
	 */
	ValueSuperbundle &reset()
	{
		for (int i = 0; i < n; ++i) memset(bundles[i].data(), 0, bundles[i].size() * sizeof(T));
		return *this;
	}

	/**
	 * @brief Fused multiply-add on all ValueBundles. 
	 * 
	 * @param rhs1 Scalar operand. 
	 * @param rhs2 Array operand. 
	 * @return ValueSuperbundle& Reference to self. 
	 */
	ValueSuperbundle &multAdd(const T &rhs1, const ValueSuperbundle<T, n> &rhs2)
	{
		for (int i = 0; i < n; ++i) bundles[i].multAdd(rhs1, rhs2.bundles[i]);
		return *this;
	}

	/**
	 * @brief Multiplication assignent. 
	 * 
	 * @param rhs Right hand side operator. 
	 * @return ValueSuperbundle& Reference to self. 
	 */
	ValueSuperbundle &operator*=(const T &rhs)
	{
		for (int i = 0; i < n; ++i) bundles[i] *= rhs;
		return *this;
	}

	/**
	 * @brief Division assignment. 
	 * 
	 * @param rhs Right hand side operator. 
	 * @return ValueSuperbundle& Reference to self. 
	 */
	ValueSuperbundle &operator/=(const T &rhs)
	{
		for (int i = 0; i < n; ++i) bundles[i] /= rhs;
		return *this;
	}

	/**
	 * @brief Addition assignment. 
	 * 
	 * @param rhs Right hand side operator. 
	 * @return ValueSuperbundle& Reference to self. 
	 */
	ValueSuperbundle &operator+=(const ValueSuperbundle<T, n> &rhs)
	{
		for (int i = 0; i < n; ++i) bundles[i] += rhs.bundles[i];
		return *this;
	}

private:
	ValueBundle<T> bundles[n]; ///< List of member ValueBundles. 
	bool hasOwnership; ///< If set to true, the ValueSuperbundle has membership of all SuperBundles' memory and deletes it upon destruction. 
};
