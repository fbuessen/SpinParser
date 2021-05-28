/**
 * @file SpinParser.hpp
 * @author Finn Lasse Buessen
 * @brief Pf-FRG solver implementation. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include "lib/Log.hpp"
#include "lib/Timestamp.hpp"
#include "lib/LoadManager.hpp"
#include "lib/Exception.hpp"
#include "FrgCommon.hpp"
#include "CommandLineOptions.hpp"
#include "TaskFileParser.hpp"

class FrgCore;

/**
 * @brief Computation status descriptor.
 */
struct ComputationStatus
{
	/**
	 * @brief Computation status identifier.
	 */
	enum struct Identifier : int
	{
		New = 0, ///< The current task is a new calculation. 
		Running = 1, ///< Calculation of task is in progress, and a checkpoint has been written. 
		Postprocessing = 2, ///< Task has finished, but measurements remain to be performed. 
		Finished = 3 ///< Task has finished. 
	};
	
	Identifier statusIdentifier; ///< Computation status identfier;
	Timestamp::Time startTime; ///< Calculation start time. 
	Timestamp::Time checkpointTime; ///< Calculation last checkpoint time. 
	Timestamp::Time endTime; ///< Calculation end time. 
};

struct Fileset
{
	std::string taskFile; ///< Path to the task file. 
	std::string obsFile; ///< Path to the observable file. 
	std::string dataFile; ///< Path to the data file used for deferred measurements. 
	std::string checkpointFile; ///< Path to the checkpoint file. 
};

/**
 * @brief Principal object and interface for the solution of pf-FRG flow equations. 
 * @details The SpinParser object provides the central interface for the solution of pf-FRG flow equations. 
 * It is designed to be a singleton, which is automatically created and retrieved via the SpinParser::spinParser() function. 
 * The solution of flow equations is then launched via the SpinParser::run() function, which is provided with the launch parameters argc and argv.
 * Once the run() function has been triggered, the SpinParser will read in the task file, prepare the specified lattice spin model, and launch the numerics core for the solution of the respective flow equations.  
 */
class SpinParser
{
public:
	/**
	 * @brief Retrieve the SpinParser singleton. 
	 * 
	 * @return SpinParser* Pointer to the singleton object. 
	 */
	static SpinParser *spinParser();

	/**
	 * @brief Launch the SpinParser. 
	 * 
	 * @param argc Launch parameter argc, as provided by the operating system. 
	 * @param argv Launch parameter argv, as provided by the operating system. 
	 * @return int Returns 0 on success, and 1 if an error occured. 
	 */
	int run(int argc, char **argv);

	/**
	 * @brief Query wheter the current instance is the MPI master rank.
	 *
	 * @return bool Returns true, if the current instance is the MPI master rank. Otherwise, returns false. 
	 */
	bool isMasterRank() const;

	/**
	 * @brief Get the current computation status. 
	 *
	 * @return ComputationStatus Computation status descriptor.
	 */
	ComputationStatus getComputationStatus() const;

	/**
	 * @brief Get file names of output files. 
	 */
	Fileset getFileset() const;

	/**
	 * @brief Retrieve the internal command line parser. 
	 * 
	 * @return CommandLineOptions* Internal command line parser. 
	 */
	CommandLineOptions *getCommandLineOptions() const;

	/**
	 * @brief Retrieve the internal task file parser. 
	 * 
	 * @return TaskFileParser* Internal task file parser. 
	 */
	TaskFileParser *getTaskFileParser() const;

	/**
	 * @brief Retrieve the internal load manager. 
	 * 
	 * @return HMP::LoadManager* Internal load manager.
	 */
	HMP::LoadManager *getLoadManager() const;

	/**
	 * @brief Retrieve the internal numerics core.
	 * 
	 * @return FrgCore* Internal numerics core.
	 */
	FrgCore *getFrgCore() const;

protected:
	/**
	 * @brief Construct a new SpinParser object. 
	 */
	SpinParser();

	/**
	 * @brief Destroy the SpinParser object. 
	 */
	~SpinParser();

	/**
	 * @brief Run the numerics core and apply the differential equation solver. 
	 */
	void runCore();

	/**
	 * @brief Write current state to checkpoint file. 
	 */
	void writeCheckpoint();

	static SpinParser *_spinParserInstance; ///< Singleton instance of the SpinParser. 
	bool _isMasterRank; ///< True, if the current instance is the MPI master rank, false otherwise. 
	ComputationStatus _computationStatus; ///< Computation status. 
	Fileset _fileset; ///< Output file names. 

	CommandLineOptions *_commandLineOptions; ///< Internal command line parser. 
	TaskFileParser *_taskFileParser; ///< Internal task file parser. 
	HMP::LoadManager *_loadManager; ///< Internal load manager. 
	FrgCore *_frgCore; ///< Internal numerics core. 
};