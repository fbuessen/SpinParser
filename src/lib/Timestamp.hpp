/**
 * @file Timestamp.hpp
 * @author Finn Lasse Buessen
 * @brief Provide formatted timestamp strings. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <boost/date_time.hpp>

namespace Timestamp
{
	/**
	 * @brief Data structure to jold time and date information. 
	 */
	typedef boost::posix_time::ptime Time;
	typedef boost::posix_time::time_duration Duration;

	/**
	 * @brief Get the current time. 
	 * 
	 * @return Time Current time. 
	 */
	inline Time time()
	{
		return boost::posix_time::second_clock::local_time();
	}

	/**
	 * @brief Get the current time plus some offset into the furutre. 
	 * 
	 * @param offset Offset duration in seconds. 
	 * @return Time Current time plus offset. 
	 */
	inline Time time(const int offset)
	{
		return boost::posix_time::second_clock::local_time() + boost::posix_time::time_duration(0, 0, offset, 0);
	}

	/**
	 * @brief Get a time object from a formatted timestamp string. 
	 * 
	 * @param timestamp Formatted timestamp string. 
	 * @return Time Time as described by the timestamp. 
	 */
	inline Time time(const std::string &timestamp)
	{
		return boost::posix_time::time_from_string(timestamp);
	}

	/**
	 * @brief Check whether the time object is older than a given number of seconds.
	 *
		@param time Time object.
		@param offset Age threshold in seconds.
		@return bool Return true if the time object is older than the specified amount of seconds. Otherwise, return false.
	 */
	inline bool isOlder(const Time &time, const int offset)
	{
		return (Timestamp::time() - time).total_seconds() > offset;
	}

	/**
	 * @brief Retrieve a timestamp string for the specified time. 
	 * 
	 * @param time Specified time. 
	 * @return std::string Timestamp string. 
	 */
	inline std::string timestamp(const Time &time)
	{
		return boost::posix_time::to_simple_string(time);
	}

	/**
	 * @brief Retrieve a timestamp string for the current time.
	 *
	 * @return std::string Timestamp string.
	 */
	inline std::string timestamp()
	{
		return timestamp(time());
	}
}
