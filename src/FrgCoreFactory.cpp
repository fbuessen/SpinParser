/**
 * @file FrgCoreFactory.cpp
 * @author Finn Lasse Buessen
 * @brief Factory function to create a new FrgCore. 
 * @details For the implementation of custom FrgCores, the FrgCoreFactory::newFrgCore routine should be extended such that it can recignize an identifier for the custom FrgCore and create an instance thereof accordingly. 
 * 
 * @copyright Copyright (c) 2020
 */

#include "lib/Exception.hpp"
#include "lib/InputParser.hpp"
#include "FrgCoreFactory.hpp"

//SU2
#include "SU2/SU2FrgCore.hpp"
#include "SU2/SU2MeasurementCorrelation.hpp"
//XYZ
#include "XYZ/XYZFrgCore.hpp"
#include "XYZ/XYZMeasurementCorrelation.hpp"
//TRI
#include "TRI/TRIFrgCore.hpp"
#include "TRI/TRIMeasurementCorrelation.hpp"


FrgCore *FrgCoreFactory::newFrgCore(const std::string &identifier, const SpinModel &model, const std::vector<MeasurementSpecification> &measurements, const std::map<std::string, std::string> &options)
{
	//create measurement objects
	std::vector<Measurement *> measurementObjects;

	for (auto specification : measurements)
	{
		//correlation measurement
		if (specification.identifier == "correlation")
		{
			Measurement *m = nullptr;
			if (identifier == "SU2") m = new SU2MeasurementCorrelation(specification.output, specification.minCutoff, specification.maxCutoff, specification.defer);
			else if (identifier == "XYZ") m = new XYZMeasurementCorrelation(specification.output, specification.minCutoff, specification.maxCutoff, specification.defer);
			else if (identifier == "TRI") m = new TRIMeasurementCorrelation(specification.output, specification.minCutoff, specification.maxCutoff, specification.defer);
			else throw Exception(Exception::Type::InitializationError, "Measurement [correlation]: Unknown model symmetry '" + identifier + "'.");

			Log::log << Log::LogLevel::Info << "Added measurement [correlation]." << Log::endl;
			measurementObjects.push_back(m);
		}
		else throw Exception(Exception::Type::InitializationError, "Measurement: Unknown measurement type '" + identifier + "'.");

	}

	if (identifier == "SU2") return new SU2FrgCore(model, measurementObjects, options);
	else if (identifier == "XYZ") return new XYZFrgCore(model, measurementObjects, options);
	else if (identifier == "TRI") return new TRIFrgCore(model, measurementObjects, options);
	else throw Exception(Exception::Type::ArgumentError, "Spin model identifier '" + identifier + "' does not exist.");
}