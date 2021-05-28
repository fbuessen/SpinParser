/**
 * @file EffectiveAction.hpp
 * @author Finn Lasse Buessen
 * @brief Virtual implementation of a flowing effective action. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <hdf5.h>
#include "lib/Log.hpp"
#include "lib/Exception.hpp"

/**
 * @brief Virtual implementation of a flowing effective action. 
 * @details Concrete implementations should implement datastructures as desired -- e.g. two-point vertex and four-point vertex information. 
 * This object provides the virtual interface accessed by FrgCore implementations and the SpinParser. 
 */
struct EffectiveAction
{
public:
	/**
	 * @brief Construct a new EffectiveAction object.
	 */
	EffectiveAction() : cutoff(0.0f) {};

	/**
	 * @brief Virtual destructor.
	 */
	virtual ~EffectiveAction() {};

	/**
	 * @brief Write all internal data to a checkpoint at the specified file path and return the identifier of the checkpoint written. The checkpoint identifier is a non-negative integer that enumerates all the checkpoint datasets in the output file, starting at zero.
	 * 
	 * @param dataFilePath Checkpoint file path. 
	 * @param append If set to false, overwrite existing checkpoint. Otherwise, append checkpoint if no previous checkpoint at the same cutoff value exists. If a checkpoint at the same cutoff value already exists, do nothing.
	 *
	 * @return int Identifier of the checkpoint written. If the writing process has been skipped, return -1. 
	 */
	virtual int writeCheckpoint(const std::string &dataFilePath, const bool append = false) const = 0;

	/**
	 * @brief Read internal data from a checkpoint with the specified checkpoint identifier at the specified file path. 
	 * If the checkpoint identifier is set to -1, read the most recent checkpoint.
	 * 
	 * @param datafilePath Checkpoint file path. 
	 * @param checkpointId Identifier of the checkpoint to read. 
	 * @return bool Return true if a checkpoint was read successfully; otherwise return false. 
	 */
	virtual bool readCheckpoint(const std::string &datafilePath, const int checkpointId = -1) = 0;

	/**
	 * @brief Indicate whether the vertex has diverged to NaN. 
	 *
	 * @return bool Return true if the vertex has diverged, otherwise return false. 
	 */
	virtual bool isDiverged() const = 0;

	float cutoff; ///< Value of the RG cutoff. 
};