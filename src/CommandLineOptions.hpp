/**
 * @file CommandLineOptions.hpp
 * @author Finn Lasse Buessen
 * @brief Parser for command line arguments. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <string>

/**
 * @brief Parser object, which can be fed with argc/argv information and which then holds the parsed values in its member variables. 
 */
class CommandLineOptions
{
public:
	/**
	 * @brief Construct a new CommandLineOptions object and parse argc/argv information. 
	 * 
	 * @param argc Argument count as passed to the program call. 
	 * @param argv Argument vector as passed to the program call. 
	 */
	CommandLineOptions(int argc, char **argv);

	/**
	 * @brief Retrieve the '--help' flag setting. 
	 * 
	 * @return bool Return true, if the '--help' flag is set. Otherwise, return false.
	 */
	bool help() const;

	/**
	 * @brief Retrieve the '--verbose' flag setting. 
	 * 
	 * @return bool Return true, if the '--verbose' flag is set. Otherwise, return false.
	 */
	bool verbose() const;

	/**
	 * @brief Retrieve the value of the '--checkpointTime' flag. 
	 * 
	 * @return int Value of the '--checkpointTime' flag. 
	 */
	int checkpointTime() const;

	/**
	 * @brief Retrieve the '--forceRestart' flag setting. 
	 * 
	 * @return bool Return true, if the '--forceRestart' flag is set. Otherwise, return false.
	 */
	bool forceRestart() const;

	/**
	 * @brief Retrieve the '--defer' flag setting. 
	 * 
	 * @return bool Return true, if the '--defer' flag is set. Otherwise, return false.
	 */
	bool deferMeasurements() const;

	/**
	 * @brief Retrieve the '--debugLattice' flag setting. 
	 * 
	 * @return bool Return true, if the '--debugLattice' flag is set. Otherwise, return false.
	 */
	bool debugLattice() const;

	/**
	 * @brief Retrieve the value of the '--taskFile' flag.
	 * 
	 * @return std::string Value of the '--taskFile' flag.
	 */
	std::string taskFile() const;

	/**
	 * @brief Retrieve the value of the '--resourcePath' flag.
	 * 
	 * @return std::string Value of the '--resourcePath' flag.
	 */
	std::string resourcePath() const;

protected:
	bool _help; ///< Help flag '--help' is set.
	bool _verbose; ///< Verbose flag '--verbose' is set.
	int _checkpointTime; ///< Value of the '--checkpointTime' argument. 
	bool _forceRestart; ///< Force flag '--forceRestart' is set. 
	bool _deferMeasurements; ///< Defer flag '--defer' is set. 
	bool _debugLattice; ///< Lattice debug flag '--debugLattice' is set. 
	std::string _taskFile; ///< Value of the '--taskFile' argument. 
	std::string _resourcePath; ///< Value of the '--resourcePath' argument. 
};