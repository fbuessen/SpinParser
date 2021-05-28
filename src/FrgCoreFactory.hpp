/**
 * @file FrgCoreFactory.hpp
 * @author Finn Lasse Buessen
 * @brief Factory function to create a new FrgCore. 
 * @details For the implementation of custom FrgCores, the FrgCoreFactory::newFrgCore routine should be extended such that it can recignize an identifier for the custom FrgCore and create an instance thereof accordingly. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <string>
#include <vector>
#include <map>
#include "FrgCore.hpp"

namespace FrgCoreFactory
{
	/**
	 * @brief Abstract specification of a measurement protocol. 
	 */
	struct MeasurementSpecification
	{
		std::string identifier; ///< String-form identifier of the measurement type, as specified in the task file. 
		std::string output; ///< Output file for the measurement results. 
		float minCutoff; ///< Minimal cutoff value for the protocol to be invoked. 
		float maxCutoff; ///< Maximal cutoff value for the protocol to be invoked. 
		bool defer; ///< Defer flag. If set to true, the measurement will only be invoked in the postprocessing stage. 
		std::vector<std::pair<std::string, std::string>> options; ///< String-form protocol modifiers as specified in the task file. 
	};

	/**
	 * @brief Create a new FrgCore for given symmetry identifier, spin model, and measurement protocols. 
	 * 
	 * @param identifier String-form symmetry identifier, as specified in the task file. 
	 * @param model Spin model used to determine the initial conditions of the flow. 
	 * @param measurements  Measurement protocols to invoke during the execution of the core. 
	 * @param options String-form core modifiers as specified in the task file. 
	 * @return FrgCore* Pointer to the new FrgCore object. 
	 */
	FrgCore *newFrgCore(const std::string &identifier, const SpinModel &model, const std::vector<MeasurementSpecification> &measurements, const std::map<std::string, std::string> &options);
}