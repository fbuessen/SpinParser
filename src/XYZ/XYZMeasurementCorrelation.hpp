/**
 * @file XYZMeasurementCorrelation.hpp
 * @author Finn Lasse Buessen
 * @brief Correlation measurement for models with diagonal interactions.
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include "Measurement.hpp"

/**
 * @brief Correlation measurement for models with diagonal interactions.
 */
class XYZMeasurementCorrelation : public Measurement
{
public:
	/**
	 * @brief Construct a new XYZMeasurementCorrelation object to measure correlations. 
	 * @see Measurement
	 * 
	 * @param outfile Filename where to write the result file. 
	 * @param minCutoff Minimum cutoff above which to invoke the measurement protocol. 
	 * @param maxCutoff Maximum cutoff below which to invoke the measurement protocol. 
	 * @param defer If set to true, measurements are deferred to the postprocessing stage. 
	 */
	XYZMeasurementCorrelation(const std::string &outfile, const float minCutoff, const float maxCutoff, const bool defer);
	
	/**
	 * @brief Destroy the XYZMeasurementCorrelation object. 
	 */
	~XYZMeasurementCorrelation();

	/**
	 * @brief Take measurement. 
	 * @see Measurement::takeMeasurement()
	 *  
	 * @param state Effective action object to perform the measurement on. 
	 * @param isMasterTask If set to true, the function call should be responsible for writing the output file. 
	 */
	void takeMeasurement(const EffectiveAction &state, const bool isMasterTask) const override;
	
private:
	/**
	 * @brief Calculate the correlation for a linear iterator in the frequency list. 
	 * 
	 * @param iterator Frequency iterator. 
	 */
	void _calculateCorrelation(const int iterator) const;

	/**
	 * @brief Write the meta information contained in the output file.
	 * @details If the output file does not yet exist, a new one is created.
	 * If the HDF5 group `observableGroup/meta` already exists, the writing process is skipped.
	 *
	 * @param observableGroup Name of the output HDF5 group.
	 */
	void _writeOutfileHeader(const std::string &observableGroup) const;

	/**
	 * @brief Write a correlation dataset of the output file. If a dateset at the current cutoff value already exists, the writing process is skipped.
	 *
	 * @param observableGroup Name of the output HDF5 group.
	 * @param correlation Correlation data to write.
	 */
	void _writeOutfileCorrelation(const std::string &observableGroup, const float *correlation) const;

	float _currentCutoff; ///< Cutoff at which the correlations have been computed. 
	float *_correlationsDD; ///< Buffer for density correlation measurements. 
	float *_correlationsXX; ///< Buffer for Sx-Sx correlation measurements. 
	float *_correlationsYY; ///< Buffer for Sy-Sy correlation measurements. 
	float *_correlationsZZ; ///< Buffer for Sz-Sz correlation measurements. 
	int _memoryStepLattice; ///< Memory stride in the correlation buffers. 
};