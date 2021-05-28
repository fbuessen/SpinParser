/**
 * @file FrgCore.hpp
 * @author Finn Lasse Buessen
 * @brief Numerics core for pf-FRG calculations. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <vector>
#include "EffectiveAction.hpp"
#include "Measurement.hpp"
#include "SpinModel.hpp"
#include "SpinParser.hpp"

class SpinParser;

/**
 * @brief Virtual implementation of a pf-FRG numerics core. 
 * @details The FrgCore represents the central numerics unit for pf-FRG calculations. 
 * It defines the interface between the differential equation solver, measurement protocols, and the concrete implementation of flow equations. 
 * The FrgCore also defines the interface for reading and writing checkpoints. 
 *
 * Every specific set of flow equations in its different symmetry-constrained form is derived from this class. 
 * Custom flow equations are implemented by subclassing the FrgCore and providing implementations of the virtual member functions. 
 * 
 * New instances of FrgCore are created by the FrgCoreFactory::newFrgCore() routine.
 */
class FrgCore
{
	friend class SpinParser;
public:
	/**
	 * @brief Invoke all associated measurement protocols. 
	 */
	void takeMeasurements() const
	{
		if (SpinParser::spinParser()->getComputationStatus().statusIdentifier == ComputationStatus::Identifier::Postprocessing)
		{
			//perform deferred measurements
			for (auto m : _measurements)
			{
				if (_flowingFunctional->cutoff <= m->maxCutoff() && _flowingFunctional->cutoff >= m->minCutoff())
				{
					if (SpinParser::spinParser()->getCommandLineOptions()->deferMeasurements() || m->isDeferred()) m->takeMeasurement(*_flowingFunctional, SpinParser::spinParser()->isMasterRank());
				}
			}
		}
		else
		{
			//perform non-deferred measurements
			for (auto m : _measurements)
			{
				if (_flowingFunctional->cutoff <= m->maxCutoff() && _flowingFunctional->cutoff >= m->minCutoff())
				{
					if (!SpinParser::spinParser()->getCommandLineOptions()->deferMeasurements() && !m->isDeferred()) m->takeMeasurement(*_flowingFunctional, SpinParser::spinParser()->isMasterRank());
				}
			}

			//write vertex output, if deferred measurements are specified
			bool postprocessingRequired = false;
			if (SpinParser::spinParser()->getCommandLineOptions()->deferMeasurements()) postprocessingRequired = true;
			for (auto m : _measurements) if (m->isDeferred()) postprocessingRequired = true;


			if (postprocessingRequired && SpinParser::spinParser()->isMasterRank()) _flowingFunctional->writeCheckpoint(SpinParser::spinParser()->getFileset().dataFile, true);
		}
	}

	/**
	 * @brief Virtual implementation of a single RG step in the solution of the flow equations. 
	 * @details The concrete implementation of the method is expected to calculate the flow equation for the current configuration in FrgCore::flowingFunctional and populate FrgCore::flow with the results. 
	 * It is not expected to make any further modifications. 
	 * 
	 * @see FrgCore::finalizeStep()
	 */
	virtual void computeStep() = 0;

	/**
	 * @brief Virtual implementation of the finalization of a single RG step in the solution of the flow equations. 
	 * @details The concrete implementation of the method is expected to update the values of FrgCore::flowingFunctional, 
	 * based on the values of the flow FrgCore::flow and the designated new value of the frequency cutoff. 
	 * 
	 * @param newCutoff New value of the cutoff. 
	 */
	virtual void finalizeStep(float newCutoff) = 0;

	/**
	 * @brief Retrieve the flowing functional.
	 *
	 * @return EffectiveAction* Flowing functional.
	 */
	EffectiveAction *flowingFunctional() const
	{
		return _flowingFunctional;
	}

	/**
	 * @brief Retrieve the vertex flow.
	 *
	 * @return EffectiveAction* Vertex flow.
	 */
	EffectiveAction *flow() const
	{
		return _flow;
	}

	/**
	 * @brief Retrieve the list of measurements.
	 *
	 * @return std::vector<Measurement *> List of measurements.
	 */
	std::vector<Measurement *> measurements() const
	{
		return _measurements;
	}

protected:
	/**
	 * @brief Construct a new FrgCore, which takes ownership of the specified measurements.
	 * @see Measurement
	 *
	 * @param measurements List of measurement protocols to invoke during the solution of the flow equations.
	 */
	FrgCore(const std::vector<Measurement *> &measurements) : _flowingFunctional(nullptr), _flow(nullptr), _measurements(measurements) {};

	/**
	 * @brief Destroy the FrgCore object and delete any associated measurement protocols.
	 */
	virtual ~FrgCore()
	{
		while (_measurements.size() > 0)
		{
			delete _measurements.back();
			_measurements.pop_back();
		}
	}

	EffectiveAction *_flowingFunctional; ///< Representation of the current state of the effective action. 
	EffectiveAction *_flow; ///< Representation of the RG flow associated with the current state of the effective action. 
	std::vector<Measurement *> _measurements; ///< List of measurement protocols to invoke throughout the solution of the flow equations. 
};