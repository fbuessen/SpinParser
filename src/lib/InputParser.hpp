/**
 * @file InputParser.hpp
 * @author Finn Lasse Buessen
 * @brief Parse input strings to numerical values. Can resolve simple multiplications and sqrt() expressions. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <iostream>
#include <sstream>
#include "boost/regex.hpp"
#include "lib/Log.hpp"

namespace InputParser
{
	/**
	 * @brief Parse input string to double. Resolves simple multiplications, divisions and sqrt() expressions. 
	 * Input may involve terms of the form x*x, x/x, sqrt(x), and longer compositions thereof, where x is a decimal number. 
	 * 
	 * @param input Input string. 
	 * @return double Parsed numerical value. 
	 */
	inline double stringToDouble(const std::string &input)
	{
		//debug output
		boost::regex number("[-\\.\\d]+");
		bool debugFlag = !boost::regex_match(input, number);
		std::string parsedString(input);

		boost::regex squareRoot("sqrt\\(([\\.\\d]+)\\)");
		boost::regex multDiv("([\\.\\d]+)([\\*/])([\\.\\d]+)");
		boost::smatch match;
		std::stringstream replacement;
		replacement << std::fixed << std::setprecision(20);
		//parse sqrt
		while (boost::regex_search(parsedString, match, squareRoot))
		{
			replacement.str("");
			replacement << match.prefix() << sqrt(std::stod(match.str(1))) << match.suffix();
			parsedString = replacement.str();
		}
		//parse multiplication and division
		while (boost::regex_search(parsedString, match, multDiv))
		{
			replacement.str("");
			replacement << match.prefix();
			if (match.str(2) == "*") replacement << std::stod(match.str(1)) * std::stod(match.str(3));
			else replacement << std::stod(match.str(1)) / std::stod(match.str(3));
			replacement << match.suffix();
			parsedString = replacement.str();
		}

		if (debugFlag) Log::log << Log::LogLevel::Debug << "parsed input string " << input << " to " << parsedString << Log::endl;
		return std::stod(parsedString);
	}

	/**
	 * @brief Parse input string to double. Resolves simple multiplications, divisions and sqrt() expressions. 
	 * Input may involve terms of the form x*x, x/x, and sqrt(x), where x is a decimal number. 
	 * 
	 * @param input Input string. 
	 * @return float Parsed numerical value. 
	 */
	inline float stringToFloat(const std::string &input)
	{
		return float(stringToDouble(input));
	}
}