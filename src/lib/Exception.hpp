/**
 * @file Exception.hpp
 * @author Finn Lasse Buessen
 * @brief Descriptor object for exceptions. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <string>
#include <sstream>
#include <exception>

/**
 * @brief Descriptor object for exceptions. 
 */
struct Exception : public std::exception
{
public:
	/**
	* @brief Exception types.
	*/
	enum struct Type
	{
		GenericError, ///< Generic error. 
		InternalError, ///< Internal error, raised when internal inconsistencies occur. 
		BadAllocation, ///< Allocation error, raised when an allocation fails. 
		OutOfRange, ///< Out of range error, raised when out of bound memory accesses occur. 
		ArgumentError, ///< Argument error, raised when a function is invoked with an invalid argument. 
		InitializationError, ///< Initialization error, raised when initialization routines fail. 
		IOError, ///< In/out error, raised when input or output routines fail. 
		AssertionFail, ///< Assertion failure, raised when an assertion fails. 
		MpiError ///< MPI error, raised when MPI communication fails. 
	};

	/**
	 * @brief Construct a new Exception object. 
	 * 
	 * @param type Type of exception to construct. 
	 */
	Exception(const Exception::Type type) noexcept
	{
		_type = type;
		switch (type)
		{
		case Exception::Type::GenericError:
			_what = "GenericError";
			break;
		case Exception::Type::InternalError:
			_what = "InternalError";
			break;
		case Exception::Type::BadAllocation:
			_what = "BadAllocation";
			break;
		case Exception::Type::OutOfRange:
			_what = "OutOfRange";
			break;
		case Exception::Type::ArgumentError:
			_what = "ArgumentError";
			break;
		case Exception::Type::InitializationError:
			_what = "InitializationError";
			break;
		case Exception::Type::IOError:
			_what = "IOError";
			break;
		case Exception::Type::AssertionFail:
			_what = "AssertionFail";
			break;
		case Exception::Type::MpiError:
			_what = "MpiError";
			break;
		default:
			_what = "UnknownType";
			break;
		}
	}

	/**
	 * @brief Construct a new Exception object. 
	 * 
	 * @tparam T Fundamental type of the exception detail descriptor. Must be convertible to a string. 
	 * @param type Type of the exception. 
	 * @param details Detailed descriptor of the exception. 
	 */
	template <class T> Exception(const Exception::Type type, T details) noexcept : Exception(type)
	{
		std::stringstream ss;
		ss << "[" << _what << "]: " << details;
		_what = ss.str();
	}

	/**
	 * @brief Destroy the Exception object. 
	 */
	~Exception() noexcept {};

	/**
	 * @brief Retrieve string form description of the exception. 
	 * 
	 * @return const char* String form description of the exception. 
	 */
	const char* what() const noexcept
	{
		return _what.c_str();
	}

	/**
	 * @brief Retrieve exception type. 
	 * 
	 * @return Exception::Type Exception type. 
	 */
	Exception::Type type() const noexcept
	{
		return _type;
	}

protected:
	Exception::Type _type; ///< Type of the exception. 
	std::string _what; ///< String form description of the exception. 
};
