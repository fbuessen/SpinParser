/**
 * @file TaskFileParser.hpp
 * @author Finn Lasse Buessen
 * @brief Task file parser routine. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <string>
#include <set>
#include <boost/property_tree/ptree.hpp>

struct FrequencyDiscretization;
struct CutoffDiscretization;
struct Lattice;
struct ComputationStatus;
class FrgCore;

/**
 * @brief Task file parser routine. The object is bound to a single task file, which it initially reads and parses, and can thereafter write to. 
 */
class TaskFileParser
{
public:
	/**
	 * @brief Construct a new TaskFileParser object and parse the specified task file. 
	 * @details During the parsing process, the TaskFileParser will allocate and initialize (with parameters as specified in the task file) FrequencyDiscretization, CutoffDiscretization, and Lattice objects; 
	 * The pointers to those objects are stored in frequency, cutoff, and lattice, respectively. 
	 * Furthermore, it allocates the FrgCore (including the FrgCore::measurements) specified in the task file and writes it to frgCore. 
	 * The computation status associated with the task file is also returned. 
	 * 
	 * @param[in] taskFilePath Path to the task file to be parsed. 
	 * @param[out] frequency Newly generated FrequencyDiscretization. 
	 * @param[out] cutoff Newly generated CutoffDiscretization. 
	 * @param[out] lattice Newly generated Lattice. 
	 * @param[out] frgCore Newly generated FRG core. 
	 * @param[out] computationStatus Computation status of the task associated with the task file. 
	 */
	TaskFileParser(const std::string &taskFilePath, FrequencyDiscretization *&frequency, CutoffDiscretization *&cutoff, Lattice *&lattice, FrgCore *&frgCore, ComputationStatus &computationStatus);

	/**
	 * @brief Write calculation status to the task file. 
	 * 
	 * @param finalize If set to true, the task status will be set to 'finished'.
	 * @param updateCheckpoint If set to true, the last checkpoint time will be set to the current time. 
	 */
	void writeTaskFile(const ComputationStatus &computationStatus);

protected: 
	/**
	 * @brief Ensure that a given property tree contains a specific node with only known required or optional elements. Throw an Exception::Type::InitializationError if the requirements are not fulfilled.
	 * 
	 * @param tree Property tree to test. 
	 * @param node Node to test. 
	 * @param requiredChildren List of required children of the node. 
	 * @param requiredAttributes List of required attributes of the node. 
	 * @param optionalChildren List of optional children of the node. 
	 * @param optionalAttributes List of optional attributes of the node.
	 */
	void _validateProperties(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &requiredChildren, const std::set<std::string> &requiredAttributes, const std::set<std::string> &optionalChildren = {}, const std::set<std::string> &optionalAttributes = {}) const;

	/**
	 * @brief Ensure that a given property tree contains a specific node and test that the node contains a list of required children. Throw an Exception::Type::InitializationError if the requirements are not fulfilled.
	 * 
	 * @param tree Property tree to test. 
	 * @param node Node to test. 
	 * @param requiredChildren List of required children of the node. 
	 * @param treePath Property path, which is prepended to the node name in a potential exception thrown.  
	 */
	void _validateRequiredChildren(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &requiredChildren, const std::string &treePath = "") const;
	
	/**
	 * @brief Ensure that a given property tree contains a specific node and test that the node contains a list of required attributes. Throw an Exception::Type::InitializationError if the requirements are not fulfilled.
	 * 
	 * @param tree Property tree to test. 
	 * @param node Node to test. 
	 * @param requiredAttributes List of required attributes of the node. 
	 * @param treePath Property path, which is prepended to the node name in a potential exception thrown.  
	 */
	void _validateRequiredAttributes(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &requiredAttributes, const std::string &treePath = "") const;
	
	/**
	 * @brief Ensure that a given property tree contains a specific node and test that the node contains only children from a list of optional values. Throw an Exception::Type::InitializationError if the requirements are not fulfilled.
	 * 
	 * @param tree Property tree to test. 
	 * @param node Node to test. 
	 * @param optionalChildren List of optional children of the node. 
	 * @param treePath Property path, which is prepended to the node name in a potential exception thrown.  
	 */
	void _validateOptionalChildren(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &optionalChildren, const std::string &treePath = "") const;
	
	/**
	 * @brief Ensure that a given property tree contains a specific node and test that the node contains only attributes from a list of optional values. Throw an Exception::Type::InitializationError if the requirements are not fulfilled.
	 * 
	 * @param tree Property tree to test. 
	 * @param node Node to test. 
	 * @param optionalAttributes List of optional attributes of the node. 
	 * @param treePath Property path, which is prepended to the node name in a potential exception thrown.  
	 */
	void _validateOptionalAttributes(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &optionalAttributes, const std::string &treePath = "") const;

	boost::property_tree::ptree _taskFile; ///< Internal property tree representation of the parsed task file. 
};