/**
 * @file Assert.hpp
 * @author Finn Lasse Buessen
 * @brief Lightweight macro library for assertions. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include "lib/Exception.hpp"

#ifdef ENABLE_ASSERTIONS
	#define ___STR(x) #x
	#define ___STRING(x) ___STR(x)
	#define __LOCATION __FILE__ ":" ___STRING(__LINE__)
#endif

/**
 * @brief Wrapper for macro overlading.
 */
#define __RESOLVE_ASSERT_OVERLOAD(_1, _2, FUNCTIONNAME, ...) FUNCTIONNAME

/**
 * @brief Ensure that the first argument is true. Optionally provide a message as the second argument, which is printed with an error message if the assertion fails. 
 */
#define ASSERT(...) __RESOLVE_ASSERT_OVERLOAD(__VA_ARGS__, __ASSERT_OVERLOAD2, __ASSERT_OVERLOAD1)(__VA_ARGS__)

#ifdef ENABLE_ASSERTIONS
	/**
	 * @brief Ensure that a statement is true. If statement is false, throw AssertionFail exception. 
	 * 
	 * @param x Statement to assess. 
	 */
	#define __ASSERT_OVERLOAD1(x) if (!(x)) throw Exception(Exception::Type::AssertionFail, "Assertion failed at " __LOCATION ". Expression " #x " should be true.");
	
	/**
	 * @brief Ensure that a statement is true. If statement is false, throw AssertionFail exception. 
	 * 
	 * @param x Statement to assess. 
	 * @param msg Message to print with the exception. 
	 */
	#define __ASSERT_OVERLOAD2(x, msg) if (!(x)) throw Exception(Exception::Type::AssertionFail, "Assertion failed at " __LOCATION ". " msg);
#else
	/**
	 * @brief Ensure that a statement is true. If statement is false, throw AssertionFail exception. 
	 * 
	 * @param x Statement to assess. 
	 */
	#define __ASSERT_OVERLOAD1(x) 

		/**
	 * @brief Ensure that a statement is true. If statement is false, throw AssertionFail exception. 
	 * 
	 * @param x Statement to assess. 
	 * @param msg Message to print with the exception. 
	 */
	#define __ASSERT_OVERLOAD2(x, msg) 
#endif