/**
 * @file XYZFrgCore.hpp
 * @author Finn Lasse Buessen
 * @brief FrgCore implementation for models with diagonal interactions.
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include "FrgCore.hpp"

/**
 * @brief FrgCore implementation for models with diagonal interactions.
 */
class XYZFrgCore : public FrgCore
{
public:
	/**
	 * @brief Construct a new XYZFrgCore, initialize with the specified spin model and add measurements. 
	 * 
	 * @param spinModel Spin model to initialize the effective action with. 
	 * @param measurements Measurements to add. 
	 * @param options String-form list of core options as provided in the task file. 
	 */
	XYZFrgCore(const SpinModel &spinModel, const std::vector<Measurement *> &measurements, const std::map<std::string, std::string> &options);
	
	/**
	 * @brief Destroy the XYZFrgCore object. 
	 */
	~XYZFrgCore();

	/**
	 * @brief Compute flow equations. 
	 */
	void computeStep() override;

	/**
	 * @brief Finalize calculation of flow equations. 
	 * 
	 * @param newCutoff New cutoff to which to extrapolate flow. 
	 */
	void finalizeStep(const float newCutoff) override;

	float normalization; ///< Energy normalization factor. 

private:
	int dataStacks[12]; ///< References to the LoadManager::DataStack. 

	/**
	 * @brief Calculate the single-particle vertex flow for a specific linear iterator, which is expanded via XYZVertexSingleParticle::expandIterator().
	 * 
	 * @param iterator Linear iterator. 
	 */
	void _calculateVertexSingleParticle(const int iterator);

	/**
	 * @brief Calculate the two-particle vertex flow for a specific linear iterator, which is expanded via XYZVertexTwoParticle::expandIterator().
	 * 
	 * @param iterator Linear iterator. 
	 */
	void _calculateVertexTwoParticle(const int iterator);
};