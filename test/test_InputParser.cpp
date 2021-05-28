#define BOOST_TEST_MODULE "InputParserTest"
#include <boost/test/included/unit_test.hpp>
#include "lib/InputParser.hpp"


BOOST_AUTO_TEST_SUITE(InputParserTest);

BOOST_AUTO_TEST_CASE(stringToDouble)
{
	std::string s;

	BOOST_CHECK_EQUAL(InputParser::stringToDouble("1.5*3.9"), 1.5 * 3.9);
	BOOST_CHECK_EQUAL(InputParser::stringToDouble("1.5/3.9"), 1.5 / 3.9);
	BOOST_CHECK_EQUAL(InputParser::stringToDouble("-1.5*3.9"), -1.5 * 3.9);
	BOOST_CHECK_EQUAL(InputParser::stringToDouble("-1.5/3.9"), -1.5 / 3.9);
	BOOST_CHECK_EQUAL(InputParser::stringToDouble("1.5*3.9/2.1"), 1.5 * 3.9 / 2.1);
	BOOST_CHECK_EQUAL(InputParser::stringToDouble("1.5/3.9*2.1"), 1.5 / 3.9 * 2.1);
	BOOST_CHECK_EQUAL(InputParser::stringToDouble("-1.5*3.9/2.1"), -1.5 * 3.9 / 2.1);
	BOOST_CHECK_EQUAL(InputParser::stringToDouble("-1.5/3.9*2.1"), -1.5 / 3.9 * 2.1);
	BOOST_CHECK_EQUAL(InputParser::stringToDouble("-1.5*sqrt(3.9)/2.1"), -1.5 * sqrt(3.9) / 2.1);
}


BOOST_AUTO_TEST_CASE(stringToFloat)
{
	std::string s;

	BOOST_CHECK_CLOSE(InputParser::stringToFloat("1.5*3.9"), 1.5 * 3.9, 1e-4);
	BOOST_CHECK_CLOSE(InputParser::stringToFloat("1.5/3.9"), 1.5 / 3.9, 1e-4);
	BOOST_CHECK_CLOSE(InputParser::stringToFloat("-1.5*3.9"), -1.5 * 3.9, 1e-4);
	BOOST_CHECK_CLOSE(InputParser::stringToFloat("-1.5/3.9"), -1.5 / 3.9, 1e-4);
	BOOST_CHECK_CLOSE(InputParser::stringToFloat("1.5*3.9/2.1"), 1.5 * 3.9 / 2.1, 1e-4);
	BOOST_CHECK_CLOSE(InputParser::stringToFloat("1.5/3.9*2.1"), 1.5 / 3.9 * 2.1, 1e-4);
	BOOST_CHECK_CLOSE(InputParser::stringToFloat("-1.5*3.9/2.1"), -1.5 * 3.9 / 2.1, 1e-4);
	BOOST_CHECK_CLOSE(InputParser::stringToFloat("-1.5/3.9*2.1"), -1.5 / 3.9 * 2.1, 1e-4);
	BOOST_CHECK_CLOSE(InputParser::stringToFloat("-1.5*sqrt(3.9)/2.1"), -1.5 * sqrt(3.9) / 2.1, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END();