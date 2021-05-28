/**
 * @file TaskFileParser.cpp
 * @author Finn Lasse Buessen
 * @brief Task file parser routine. 
 * 
 * @copyright Copyright (c) 2020
 */

#include "TaskFileParser.hpp"
#include <set>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "lib/InputParser.hpp"
#include "lib/Timestamp.hpp"
#include "LatticeModelFactory.hpp"
#include "FrgCoreFactory.hpp"
#include "SpinParser.hpp"


TaskFileParser::TaskFileParser(const std::string &taskFilePath, FrequencyDiscretization *&frequency, CutoffDiscretization *&cutoff, Lattice *&lattice, FrgCore *&frgCore, ComputationStatus &computationStatus)
{
	//parse xml document
	boost::property_tree::read_xml(taskFilePath, _taskFile, boost::property_tree::xml_parser::no_concat_text);

	//validate global task file structure
	_validateProperties(_taskFile, "", { "task" }, {});
	_validateProperties(_taskFile, "task", { "parameters" }, {}, { "measurements", "calculation" });
	_validateProperties(_taskFile, "task.parameters", { "frequency", "cutoff", "lattice", "model" }, {});

	//computation status
	#pragma region computation status
	if (SpinParser::spinParser()->getCommandLineOptions()->forceRestart()) computationStatus.statusIdentifier = ComputationStatus::Identifier::New;
	else
	{
		if (_taskFile.get_optional<std::string>("task.calculation.<xmlattr>.status"))
		{
			std::string status = _taskFile.get<std::string>("task.calculation.<xmlattr>.status");

			if (status == "running")
			{
				computationStatus.statusIdentifier = ComputationStatus::Identifier::Running;
				
				if (_taskFile.get_optional<std::string>("task.calculation.<xmlattr>.startTime")) computationStatus.startTime = Timestamp::time(_taskFile.get<std::string>("task.calculation.<xmlattr>.startTime"));
				else throw Exception(Exception::Type::InitializationError, "Checkpoint is corrupted [start time is not specified].");
			}
			else if (status == "postprocessing")
			{
				computationStatus.statusIdentifier = ComputationStatus::Identifier::Postprocessing;
				
				if (_taskFile.get_optional<std::string>("task.calculation.<xmlattr>.startTime")) computationStatus.startTime = Timestamp::time(_taskFile.get<std::string>("task.calculation.<xmlattr>.startTime"));
				else throw Exception(Exception::Type::InitializationError, "Checkpoint is corrupted [start time is not specified].");

				if (_taskFile.get_optional<std::string>("task.calculation.<xmlattr>.checkpointTime")) computationStatus.checkpointTime = Timestamp::time(_taskFile.get<std::string>("task.calculation.<xmlattr>.checkpointTime"));
				else throw Exception(Exception::Type::InitializationError, "Checkpoint is corrupted [checkpoint time is not specified].");
			}
			else if (status == "finished")
			{
				computationStatus.statusIdentifier = ComputationStatus::Identifier::Finished;

				if (_taskFile.get_optional<std::string>("task.calculation.<xmlattr>.startTime")) computationStatus.startTime = Timestamp::time(_taskFile.get<std::string>("task.calculation.<xmlattr>.startTime"));
				else throw Exception(Exception::Type::InitializationError, "Checkpoint is corrupted [start time is not specified].");

				if (_taskFile.get_optional<std::string>("task.calculation.<xmlattr>.checkpointTime")) computationStatus.checkpointTime = Timestamp::time(_taskFile.get<std::string>("task.calculation.<xmlattr>.checkpointTime"));
				else throw Exception(Exception::Type::InitializationError, "Checkpoint is corrupted [checkpoint time is not specified].");

				if (_taskFile.get_optional<std::string>("task.calculation.<xmlattr>.endTime")) computationStatus.endTime = Timestamp::time(_taskFile.get<std::string>("task.calculation.<xmlattr>.endTime"));
				else throw Exception(Exception::Type::InitializationError, "Cannot continue calculation of task. Checkpoint seems corrupted (end time is not specified).");
			}
			else throw Exception(Exception::Type::InitializationError, "Checkpoint is corrupted [end time is not specified].");
		}
		else computationStatus.statusIdentifier = ComputationStatus::Identifier::New;
	}

	if (computationStatus.statusIdentifier == ComputationStatus::Identifier::New)
	{
		//make sure no previous checkpoint or data files exist
		if (SpinParser::spinParser()->isMasterRank())
		{
			if (boost::filesystem::exists(SpinParser::spinParser()->getFileset().checkpointFile))
			{
				if (SpinParser::spinParser()->getCommandLineOptions()->forceRestart())
				{
					Log::log << Log::LogLevel::Warning << "Checkpoint file [" + SpinParser::spinParser()->getFileset().checkpointFile + "] already exists. File will be overwritten." << Log::endl;
					boost::filesystem::remove(SpinParser::spinParser()->getFileset().checkpointFile);
				}
				else throw Exception(Exception::Type::IOError, "Checkpoint file [" + SpinParser::spinParser()->getFileset().checkpointFile + "] already exists. Specify the --forceRestart option to overwrite existing files.");
			}

			if (boost::filesystem::exists(SpinParser::spinParser()->getFileset().dataFile))
			{
				if (SpinParser::spinParser()->getCommandLineOptions()->forceRestart())
				{
					Log::log << Log::LogLevel::Warning << "Data file [" + SpinParser::spinParser()->getFileset().dataFile + "] already exists. File will be overwritten." << Log::endl;
					boost::filesystem::remove(SpinParser::spinParser()->getFileset().dataFile);
				}
				else throw Exception(Exception::Type::IOError, "Data file [" + SpinParser::spinParser()->getFileset().dataFile + "] already exists. Specify the --forceRestart option to overwrite existing files.");
			}
		}
	}
	#pragma endregion

	//frequency
	#pragma region frequency
	_validateProperties(_taskFile, "task.parameters.frequency", {}, { "discretization" }, { "min", "max", "count", "value" });

	if (_taskFile.get<std::string>("task.parameters.frequency.<xmlattr>.discretization") == "exponential")
	{
		_validateProperties(_taskFile, "task.parameters.frequency", { "min", "max", "count" }, { "discretization" });

		//populate discretization automatically
		float min = InputParser::stringToFloat(_taskFile.get<std::string>("task.parameters.frequency.min.<xmltext>"));
		if (min <= 0) throw Exception(Exception::Type::InitializationError, "Invalid task file. Parameter 'task.parameters.frequency.min' must be positive");

		float max = InputParser::stringToFloat(_taskFile.get<std::string>("task.parameters.frequency.max.<xmltext>"));
		if (max <= min) throw Exception(Exception::Type::InitializationError, "Invalid task file. Parameter 'task.parameters.frequency.max' must be greater than 'task.parameters.frequency.min'");

		float count = InputParser::stringToFloat(_taskFile.get<std::string>("task.parameters.frequency.count.<xmltext>"));
		if (count <= 1) throw Exception(Exception::Type::InitializationError, "Invalid task file. Parameter 'task.parameters.frequency.count' must be greater than 2.");

		std::vector<float> frequencies;
		float step = powf(max / min, 1.0f/(count - 1));
		for (int i = 0; i < count; ++i) frequencies.push_back(min*powf(step,float(i)));

		frequency = new FrequencyDiscretization(frequencies);
		Log::log << Log::LogLevel::Info << "Generated exponential frequency discretization with " << frequencies.size() << " values" << Log::endl;
	}
	else if (_taskFile.get<std::string>("task.parameters.frequency.<xmlattr>.discretization") == "manual")
	{
		_validateProperties(_taskFile, "task.parameters.frequency", {}, { "discretization" }, { "value" });

		std::vector<float> frequencies;
		for (auto node : _taskFile.get_child("task.parameters.frequency"))
		{
			if (node.first == "value")
			{
				if (!node.second.get_optional<std::string>("<xmltext>")) throw Exception(Exception::Type::InitializationError, "Invalid task file. Unspecified parameter value (task.parameters.frequency.value)");
				frequencies.push_back(InputParser::stringToFloat(node.second.get<std::string>("<xmltext>")));
			}
		}
		std::sort(frequencies.begin(), frequencies.end());
		for (unsigned int i = 0; i < frequencies.size() - 1; ++i)
		{
			if (frequencies[i] == frequencies[i + 1]) throw Exception(Exception::Type::InitializationError, "Invalid task file. Duplicate value in task.parameters.frequency.value");
		}

		//set up frequency parametrization
		frequency = new FrequencyDiscretization(frequencies);
	}
	else throw Exception(Exception::Type::InitializationError, "Invalid task file. Unknown attribute value '" + _taskFile.get<std::string>("task.parameters.frequency.<xmlattr>.discretization") + "' (task.parameters.frequency.discretization)");
	#pragma endregion

	//cutoff
	#pragma region cutoff
	_validateProperties(_taskFile, "task.parameters.cutoff", {}, { "discretization" }, { "min", "max", "step", "value" });

	if (_taskFile.get<std::string>("task.parameters.cutoff.<xmlattr>.discretization") == "exponential")
	{
		_validateProperties(_taskFile, "task.parameters.cutoff", { "min", "max", "step" }, { "discretization" });

		//populate discretization automatically
		float min = InputParser::stringToFloat(_taskFile.get<std::string>("task.parameters.cutoff.min.<xmltext>"));
		if (min <= 0) throw Exception(Exception::Type::InitializationError, "Invalid task file. Parameter 'task.parameters.cutoff.min' must be positive");

		float max = InputParser::stringToFloat(_taskFile.get<std::string>("task.parameters.cutoff.max.<xmltext>"));
		if (max <= min) throw Exception(Exception::Type::InitializationError, "Invalid task file. Parameter 'task.parameters.cutoff.max' must be greater than 'task.parameters.cutoff.min'");

		float step = InputParser::stringToFloat(_taskFile.get<std::string>("task.parameters.cutoff.step.<xmltext>"));
		if (step <= 0 || step >= 1) throw Exception(Exception::Type::InitializationError, "Invalid task file. Parameter 'task.parameters.cutoff.step' must be in the range (0,1)");

		std::vector<float> cutoffValues;
		while (max > min)
		{
			cutoffValues.push_back(max);
			max *= step;
		}

		cutoff = new CutoffDiscretization(cutoffValues);
		Log::log << Log::LogLevel::Info << "Generated exponential cutoff discretization with " << cutoffValues.size() << " values" << Log::endl;
	}
	else if (_taskFile.get<std::string>("task.parameters.cutoff.<xmlattr>.discretization") == "manual")
	{
		_validateProperties(_taskFile, "task.parameters.cutoff", {}, { "discretization" }, { "value" });

		//populate discretization manually
		std::vector<float> cutoffValues;
		for (auto node : _taskFile.get_child("task.parameters.cutoff"))
		{
			if (node.first == "value")
			{
				if (!node.second.get_optional<std::string>("<xmltext>")) throw Exception(Exception::Type::InitializationError, "Invalid task file. Unspecified parameter value (task.parameters.cutoff.value)");
				cutoffValues.push_back(InputParser::stringToFloat(node.second.get<std::string>("<xmltext>")));
			}
		}
		std::sort(cutoffValues.begin(), cutoffValues.end(), std::greater<float>());
		cutoff = new CutoffDiscretization(cutoffValues);
	}
	else throw Exception(Exception::Type::InitializationError, "Invalid task file. Unknown attribute value '" + _taskFile.get<std::string>("task.parameters.cutoff.<xmlattr>.discretization") + "' (task.parameters.cutoff.discretization)");
	#pragma endregion

	//lattice model
	#pragma region lattice model
	_validateProperties(_taskFile, "task.parameters.lattice", {}, { "name", "range" });
	_validateRequiredAttributes(_taskFile, "task.parameters.model", { "name", "symmetry" });

	int latticerange = std::stoi(_taskFile.get<std::string>("task.parameters.lattice.<xmlattr>.range"));
	std::string latticeName = _taskFile.get<std::string>("task.parameters.lattice.<xmlattr>.name");
	std::string modelName = _taskFile.get<std::string>("task.parameters.model.<xmlattr>.name");
	std::string coreIdentifier = _taskFile.get<std::string>("task.parameters.model.<xmlattr>.symmetry");

	std::map<std::string, std::string> options;
	for (auto node : _taskFile.get_child("task.parameters.model"))
	{
		std::string optionName = node.first;
		if (optionName == "<xmlattr>" || optionName == "<xmltext>" || optionName == "<xmlcomment>") continue;
		if (!node.second.get_optional<std::string>("<xmltext>")) throw Exception(Exception::Type::InitializationError, "Invalid task file. Unspecified parameter value (task.parameters.model." + optionName + ")");
		std::string optionValue = node.second.get<std::string>("<xmltext>");
		auto status = options.insert(std::pair<std::string, std::string>(optionName, optionValue));
		if (status.second == false) throw Exception(Exception::Type::InitializationError, "Invalid task file. Multiple definitions of parameter 'task.parameters.model." + optionName + "'");
	}

	LatticeModelFactory::LatticeUnitCell factoryLatticeUC(latticeName, SpinParser::spinParser()->getCommandLineOptions()->resourcePath());
	LatticeModelFactory::SpinModelUnitCell factorySpinUC(modelName, SpinParser::spinParser()->getCommandLineOptions()->resourcePath(), options);
	std::pair<Lattice *, SpinModel *> factoryProduct = LatticeModelFactory::newLatticeModel(factoryLatticeUC, factorySpinUC, latticerange, boost::filesystem::path(taskFilePath).replace_extension("ldf").string());
	lattice = factoryProduct.first;
	SpinModel *spinModel = factoryProduct.second;

	Log::log << Log::LogLevel::Info << Log::LogLevel::Info << "Generated lattice model." << Log::endl;
	#pragma endregion
	
	//measurements
	#pragma region measurements
	std::vector<FrgCoreFactory::MeasurementSpecification> measurements;
	if (_taskFile.get_child_optional("task.measurements"))
	{
		_validateProperties(_taskFile, "task.measurements", {}, {}, { "measurement" });

		for (auto node : _taskFile.get_child("task.measurements"))
		{
			if (node.first != "measurement") continue;
			auto measurementTask = node.second;

			_validateRequiredAttributes(measurementTask, "", { "name" }, "task.measurements.measurement");
			_validateOptionalAttributes(measurementTask, "", { "name", "output", "method", "minCutoff", "maxCutoff" }, "task.measurements.measurement");

			//observable name
			std::string measurementIdentifier = measurementTask.get<std::string>("<xmlattr>.name");

			//observable file path
			std::string output;
			if (measurementTask.get_optional<std::string>("<xmlattr>.output")) output = boost::filesystem::path(taskFilePath).remove_filename().append(measurementTask.get<std::string>("<xmlattr>.output")).string();
			else output = SpinParser::spinParser()->getFileset().obsFile;
			
			if (computationStatus.statusIdentifier == ComputationStatus::Identifier::New)
			{
				if (SpinParser::spinParser()->isMasterRank())
				{
					if (boost::filesystem::exists(output))
					{
						if (SpinParser::spinParser()->getCommandLineOptions()->forceRestart())
						{
							Log::log << Log::LogLevel::Warning << "Observable file [" << output << "] already exists. File will be overwritten. " << Log::endl;
							boost::filesystem::remove(output);
						}
						else throw Exception(Exception::Type::IOError, "Observable file [" + output + "] already exists. Specify the --forceRestart option to overwrite existing files.");
					}

				}
			}

			//deferred measurement
			bool defer = false;
			if (measurementTask.get_optional<std::string>("<xmlattr>.method"))
			{
				if (measurementTask.get<std::string>("<xmlattr>.method") == "defer") defer = true;
			}

			//min cutoff
			float minCutoff = 0.0f;
			if (measurementTask.get_optional<std::string>("<xmlattr>.minCutoff")) minCutoff = InputParser::stringToFloat(measurementTask.get<std::string>("<xmlattr>.minCutoff"));

			//max cutoff
			float maxCutoff = INFINITY;
			if (measurementTask.get_optional<std::string>("<xmlattr>.maxCutoff")) maxCutoff = InputParser::stringToFloat(measurementTask.get<std::string>("<xmlattr>.maxCutoff"));

			//read additional options
			std::vector<std::pair<std::string, std::string>> measurementOptions;
			for (auto node : measurementTask)
			{
				std::string optionName = node.first;
				if (optionName == "<xmlattr>" || optionName == "<xmltext>" || optionName == "<xmlcomment>") continue;
				if (!node.second.get_optional<std::string>("<xmltext>")) throw Exception(Exception::Type::InitializationError, "Invalid task file. Unspecified parameter value (task.measurements.measurement." + optionName + ")");
				std::string optionValue = node.second.get<std::string>("<xmltext>");
				measurementOptions.push_back(std::pair<std::string, std::string>(optionName, optionValue));
			}

			//make measurement
			FrgCoreFactory::MeasurementSpecification s;
			s.identifier = measurementIdentifier;
			s.output = output;
			s.minCutoff = minCutoff;
			s.maxCutoff = maxCutoff;
			s.defer = defer;
			s.options = measurementOptions;
			measurements.push_back(s);
		}
	}
	#pragma endregion

	//FRG core
	#pragma region FRG core
	//make frg core
	std::map<std::string, std::string> coreOptions;
	for (auto option : options)
	{
		if (std::find(spinModel->interactionParameters.begin(), spinModel->interactionParameters.end(), option.first) == spinModel->interactionParameters.end()) coreOptions.insert(option);
	}

	frgCore = FrgCoreFactory::newFrgCore(coreIdentifier, *spinModel, measurements, coreOptions);

	Log::log << Log::LogLevel::Info << Log::LogLevel::Info << "Generated FRG core with identifier " << coreIdentifier << "." << Log::endl;
	#pragma endregion

	//debugLattice output
	#pragma region debugLattice output
	if (SpinParser::spinParser()->getCommandLineOptions()->debugLattice())
	{
		//print general information
		Log::log << Log::LogLevel::Info << "Printing lattice debug information." << Log::endl;
		
		Log::log << Log::LogLevel::Info << "\tNumber of basis sites per unit cell: " << FrgCommon::lattice()._basis.size() << Log::endl;
		Log::log << Log::LogLevel::Info << "\tNumber of allocated sites: " << FrgCommon::lattice().end() - FrgCommon::lattice().begin() << " (" << (FrgCommon::lattice().end() - FrgCommon::lattice().begin()) * sizeof(int) << " bytes)" << Log::endl;
		int inRangeCount = 0;
		for (auto i1 = FrgCommon::lattice().getRange(0); i1 != FrgCommon::lattice().end(); ++i1) ++inRangeCount;
		Log::log << Log::LogLevel::Info << "\tNumber of sites within range of origin: " << inRangeCount << Log::endl;
		Log::log << Log::LogLevel::Info << "\tNumber of parametrized sites: " << FrgCommon::lattice().size << Log::endl;

		Log::log << Log::LogLevel::Info << "\tParametrization is as follows: " << Log::endl;
		Log::log << Log::LogLevel::Info << "\t\tid: (a0, a1, a2, b)" << Log::endl;
		for (int i = 0; i < FrgCommon::lattice().size; ++i)
		{
			auto p = FrgCommon::lattice().getSiteParameters(FrgCommon::lattice().fromParametrization(i));
			Log::log << Log::LogLevel::Info << "\t\t" << i << ": (" << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << ", " << std::get<3>(p) << ")" << Log::endl;
		}
	}
	#pragma endregion
}

void TaskFileParser::writeTaskFile(const ComputationStatus &computationStatus)
{
	//remove previous calculation log from taskfile
	if (_taskFile.get_child_optional("task.calculation")) _taskFile.get_child("task").erase("calculation");
	
	//generate new calculation log
	if (computationStatus.statusIdentifier != ComputationStatus::Identifier::New)
	{
		_taskFile.put("task.calculation.<xmlattr>.startTime", Timestamp::timestamp(computationStatus.startTime));
		_taskFile.put("task.calculation.<xmlattr>.checkpointTime", Timestamp::timestamp(computationStatus.checkpointTime));

		if (computationStatus.statusIdentifier == ComputationStatus::Identifier::Running) _taskFile.put("task.calculation.<xmlattr>.status", "running");
		else if (computationStatus.statusIdentifier == ComputationStatus::Identifier::Postprocessing) _taskFile.put("task.calculation.<xmlattr>.status", "postprocessing");
		else if (computationStatus.statusIdentifier == ComputationStatus::Identifier::Finished)
		{
			_taskFile.put("task.calculation.<xmlattr>.endTime", Timestamp::timestamp(computationStatus.endTime));
			_taskFile.put("task.calculation.<xmlattr>.status", "finished");
		}
	}
	
	//write file
	std::ofstream file(SpinParser::spinParser()->getFileset().taskFile, std::ios::out);
	if (!file.is_open()) throw Exception(Exception::Type::IOError, "Could not open task file for writing");
	boost::property_tree::write_xml(file, _taskFile);
		file.close();
}

void TaskFileParser::_validateProperties(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &requiredChildren, const std::set<std::string> &requiredAttributes, const std::set<std::string> &optionalChildren, const std::set<std::string> &optionalAttributes) const
{
	_validateRequiredChildren(tree, node, requiredChildren);
	_validateRequiredAttributes(tree, node, requiredAttributes);

	std::set<std::string> allChildren(requiredChildren);
	allChildren.insert(optionalChildren.begin(), optionalChildren.end());
	_validateOptionalChildren(tree, node, allChildren);
	std::set<std::string> allAttributes(requiredAttributes);
	allAttributes.insert(optionalAttributes.begin(), optionalAttributes.end());
	_validateOptionalAttributes(tree, node, allAttributes);
}

void TaskFileParser::_validateRequiredChildren(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &requiredChildren, const std::string &treePath) const
{
	if (!tree.get_child_optional(node)) throw Exception(Exception::Type::InitializationError, "Invalid task file. Missing parameter '" + node + "'.");

	//check if all required children exist
	for (auto v : requiredChildren) if (!tree.get_child_optional(((node == "") ? "" : node + ".") + v)) throw Exception(Exception::Type::InitializationError, "Invalid task file. Missing parameter '" + ((treePath == "") ? "" : treePath + ".") + ((node == "") ? "" : node + ".") + v + "'.");
}

void TaskFileParser::_validateRequiredAttributes(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &requiredAttributes, const std::string &treePath) const
{
	if (!tree.get_child_optional(node)) throw Exception(Exception::Type::InitializationError, "Invalid task file. Missing parameter '" + node + "'.");

	//check if all required attributes exist
	for (auto v : requiredAttributes) if (!tree.get_child_optional(((node == "") ? "" : node + ".") + "<xmlattr>." + v)) throw Exception(Exception::Type::InitializationError, "Invalid task file. Missing attribute '" + ((treePath == "") ? "" : treePath + ".") + ((node == "") ? "" : node + ".") + v + "'.");
}

void TaskFileParser::_validateOptionalChildren(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &optionalChildren, const std::string &treePath) const
{
	if (!tree.get_child_optional(node)) throw Exception(Exception::Type::InitializationError, "Invalid task file. Missing parameter '" + node + "'.");

	//check if any unknown child appears
	std::set<std::string> allChildren(optionalChildren);
	allChildren.insert({ "<xmlattr>", "<xmlcomment>", "<xmltext>" });
	for (auto v : tree.get_child(node)) if (allChildren.find(v.first) == allChildren.end()) throw Exception(Exception::Type::InitializationError, "Invalid task file. Unknown parameter '" + ((treePath == "") ? "" : treePath + ".") + ((node == "") ? "" : node + ".") + v.first + "'.");
}

void TaskFileParser::_validateOptionalAttributes(const boost::property_tree::ptree &tree, const std::string &node, const std::set<std::string> &optionalAttributes, const std::string &treePath) const
{
	if (!tree.get_child_optional(node)) throw Exception(Exception::Type::InitializationError, "Invalid task file. Missing parameter '" + node + "'.");

	//check if any unknown attribute appears
	if (tree.get_child_optional(((node == "") ? "" : node + ".") + "<xmlattr>"))
	{
		for (auto v : tree.get_child(((node == "") ? "" : node + ".") + "<xmlattr>")) if (optionalAttributes.find(v.first) == optionalAttributes.end()) throw Exception(Exception::Type::InitializationError, "Invalid task file. Unknown attribute '" + ((treePath == "") ? "" : treePath + ".") + ((node == "") ? "" : node + ".") + v.first + "'.");
	}
}