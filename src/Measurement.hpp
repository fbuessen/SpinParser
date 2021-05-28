/**
 * @file Measurement.hpp
 * @author Finn Lasse Buessen
 * @brief Abstract measurement protocol definition for pf-FRG calculations. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <string>
#include <vector>
#include "lib/LoadManager.hpp"

struct EffectiveAction;

/**
 * @brief Virtual implementation of a measurement protocol. 
 * @details Measurement protocols are part of the FrgCore. Concrete implementations of the protocol may take specific measurements. 
 * A measurement may potentially have two special properties; It can be a deferred measurement, and it can be a load managed measurement. 
 * 
 * A deferred measurement will not compute any measurements during the solution of the flow equations. 
 * Instead, the measurements are copmuted in the Postprocessing phase, which is entered by re-running the SpinParser with the same task file that includes deferred measurements. 
 * In order to enable postprocessing measurements, during the solution of the flow equations, the vertex data itself is written to the disk. 
 * The `deferred` property can be chosen at initialization time of the measurement. 
 * 
 * A load managed measurement does not necessarily run all computations on its own. 
 * Instead, it can provide a list of LoadManager::DataStack ids, which are then calculated in the FrgCore::computeStep() phase. 
 * This allows the FrgCore to perform better load balancing between the different flow equations and the calculations required for the measurements. 
 * The `load managed` propery is inherent to the measurement. Derived measurement classes should initialize the member variable Meausrement::_isLoadManaged with the desired value. 
 */
class Measurement
{
public:
	/**
	 * @brief Construct a non load managed Measurement object.
	 *
	 * @param outfile Filename where to write the result file.
	 * @param minCutoff Minimum cutoff above which to invoke the measurement protocol.
	 * @param maxCutoff Maximum cutoff below which to invoke the measurement protocol.
	 * @param isDeferred If set to true, measurements are deferred to the postprocessing stage.
	 */
	Measurement(const std::string &outfile, const float minCutoff, const float maxCutoff, const bool isDeferred);

	/**
	 * @brief Virtual destructor. 
	 */
	virtual ~Measurement();
	
	/**
	 * @brief Virtual implementation of the measurement routine. This routine is called from the SpinParser whenever a measurement should be performed.
	 * Measurement sould be performed on the specified effective action.
	 *
	 * @param state Effective action object to perform the measurement on.
	 * @param isMasterTask If set to true, the function call should be responsible for writing the output file.
	 */
	virtual void takeMeasurement(const EffectiveAction &state, const bool isMasterTask) const = 0;

	/**
	 * @brief Return the filename of the output file.
	 *
	 * @return std::string Output filename.
	 */
	std::string outfile() const;

	/**
	 * @brief Return the minimum cutoff above which the measurement protocol is invoked.
	 *
	 * @return float Minimum cutoff value.
	 */
	float minCutoff() const;

	/**
	 * @brief Return the maximum cutoff value below which the measurement protocol is invoked.
	 *
	 * @return float Maximum cutoff value.
	 */
	float maxCutoff() const;

	/**
	 * @brief Query whether the measurement protocol is a deferred measurement.
	 *
	 * @return bool Return true, if the measurement is a deferred measurement. Return false otherwise.
	 */
	bool isDeferred() const;

	/**
	 * @brief Query whether the measurement protocol is load managed.
	 * @details A measurement protocol which is load managed provides a list of LoadManager::DataStack ids in the Measurement::getLoadManagedStacks() routine.
	 * Those stacks should then be calculated in the FrgCore::computeStep() function; This mechanism allows for better load balancing.
	 *
	 * @return bool Return true, if the measurement is load managed. Return false otherwise.
	 */
	bool isLoadManaged() const;

	/**
	 * @brief Return a list of LoadManager::DataStack ids to compute in the FrgCore::computeStep() function.
	 * @see Measurement::isLoadManaged()
	 *
	 * @return std::vector<int> List of DataStack ids.
	 */
	std::vector<HMP::StackIdentifier> getLoadManagedStacks() const;

protected:
	/**
	 * @brief Construct a new Measurement object.
	 *
	 * @param outfile Filename where to write the result file.
	 * @param minCutoff Minimum cutoff above which to invoke the measurement protocol.
	 * @param maxCutoff Maximum cutoff below which to invoke the measurement protocol.
	 * @param isDeferred If set to true, measurements are deferred to the postprocessing stage.
	 * @param isLoadManaged If set to true, measurement is assumed to be load managed.
	 */
	Measurement(const std::string &outfile, const float minCutoff, const float maxCutoff, const bool isDeferred, const bool isLoadManaged);

	bool _isLoadManaged; ///< If set to true, the measurement protocol is considered to be load managed. Derived classes should initialize this variable with the desired value in the constructor. 
	std::vector<HMP::StackIdentifier> _loadManagedStacks; ///< Contains a list of load managed stack identifiers. Derived classis should initialize this list in the constructor. 

private:
	std::string _outfile; ///< Filename where to write the result file.
	float _minCutoff; ///< Minimum cutoff above which to invoke the measurement protocol. 
	float _maxCutoff; ///< Maximum cutoff below which to invoke the measurement protocol. 
	bool _isDeferred; ///< If set to true, measurements are deferred to the postprocessing stage. 
};