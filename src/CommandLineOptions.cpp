/**
 * @file CommandLineOptions.cpp
 * @author Finn Lasse Buessen
 * @brief Parser for command line arguments. 
 * 
 * @copyright Copyright (c) 2020
 */

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "CommandLineOptions.hpp"

CommandLineOptions::CommandLineOptions(int argc, char **argv)
{
	namespace po = boost::program_options;
	
	po::options_description generalOptions("General options");
	generalOptions.add_options()
		("help,h", po::bool_switch(), "print help message and exit")
		("resourcePath,r", po::value<std::string>()->value_name("DIR"), "search path for .xml resource files");

	po::options_description checkpointingOptions("Checkpointing options");
	checkpointingOptions.add_options()
		("checkpointTime,t", po::value<int>()->default_value(3600)->value_name("TIME"), "checkpoint interval in seconds")
		("forceRestart,f", po::bool_switch(), "start new calculation even if checkpoint data is available")
		("defer,d", po::bool_switch(), "archive all vertex data for deferred measurements in post processing");

	po::options_description outputOptions("Output options");
	outputOptions.add_options()
		("verbose,v", po::bool_switch(), "enable verbose output")
		("debugLattice", po::bool_switch(), "print lattice debug information in .ldf format");

	po::options_description hiddenOptions("Hidden options");
	hiddenOptions.add_options()
		("taskFile", po::value<std::string>()->required(), "taskFile");

	po::options_description allOptions;
	allOptions.add(generalOptions).add(checkpointingOptions).add(outputOptions).add(hiddenOptions);

	po::positional_options_description positionalOptions;
	positionalOptions.add("taskFile", -1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(allOptions).positional(positionalOptions).run(), vm);
	
	_help = vm["help"].as<bool>();
	_verbose = vm["verbose"].as<bool>();
	_checkpointTime = vm["checkpointTime"].as<int>();
	_forceRestart = vm["forceRestart"].as<bool>();
	_deferMeasurements = vm["defer"].as<bool>();
	_debugLattice = vm["debugLattice"].as<bool>();
	_taskFile = (vm.count("taskFile")) ? vm["taskFile"].as<std::string>() : "";
	if (vm.count("resourcePath")) _resourcePath = vm["resourcePath"].as<std::string>();
	else
	{
		if (boost::filesystem::is_directory(boost::filesystem::path(argv[0]).remove_filename().parent_path().append("res").string())) _resourcePath = boost::filesystem::path(argv[0]).remove_filename().parent_path().append("res").string();
		else if (boost::filesystem::is_directory(boost::filesystem::path(argv[0]).remove_filename().parent_path().parent_path().append("res").string())) _resourcePath = boost::filesystem::path(argv[0]).remove_filename().parent_path().parent_path().append("res").string();
		else _resourcePath = boost::filesystem::path(argv[0]).remove_filename().string();
	}

	if (_help)
	{
		std::cout << "NAME" << std::endl;
		std::cout << "\tSpinParser - Spin Pseudofermion Algorithms for Research on Spin Ensembles via Renormalization" << std::endl << std::endl;
		std::cout << "SYNOPSIS" << std::endl;
		std::cout << "\tSpinParser [OPTION]... FILE" << std::endl << std::endl;
		std::cout << "DESCRIPTION" << std::endl;
		std::cout << "\tThe SpinParser allows to solve pf-FRG flow equations with model parameters specified in FILE. " << std::endl << std::endl;
		std::cout << "\tMandatory arguments to long options are mandatory for short options too. " << std::endl << std::endl;

		std::cout << generalOptions << std::endl << outputOptions << std::endl << checkpointingOptions << std::endl;
	}
	else po::notify(vm);
}

bool CommandLineOptions::help() const
{
	return _help;
}

bool CommandLineOptions::verbose() const
{
	return _verbose;
}

int CommandLineOptions::checkpointTime() const
{
	return _checkpointTime;
}

bool CommandLineOptions::forceRestart() const
{
	return _forceRestart;
}

bool CommandLineOptions::deferMeasurements() const
{
	return _deferMeasurements;
}

bool CommandLineOptions::debugLattice() const
{
	return _debugLattice;
}

std::string CommandLineOptions::taskFile() const
{
	return _taskFile;
}

std::string CommandLineOptions::resourcePath() const
{
	return _resourcePath;
}
