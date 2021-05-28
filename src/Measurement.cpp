/**
 * @file Measurement.cpp
 * @author Finn Lasse Buessen
 * @brief Abstract measurement protocol definition for pf-FRG calculations. 
 * 
 * @copyright Copyright (c) 2020
 */

#include "Measurement.hpp"
#include "lib/Exception.hpp"
#include "FrgCommon.hpp"

Measurement::Measurement(const std::string &outfile, const float minCutoff, const float maxCutoff, const bool isDeferred) : _isLoadManaged(false), _outfile(outfile), _minCutoff(minCutoff), _maxCutoff(maxCutoff), _isDeferred(isDeferred) {}

Measurement::Measurement(const std::string &outfile, const float minCutoff, const float maxCutoff, const bool isDeferred, const bool isLoadManaged) : _isLoadManaged(isLoadManaged), _outfile(outfile), _minCutoff(minCutoff), _maxCutoff(maxCutoff), _isDeferred(isDeferred) {}

Measurement::~Measurement() {};

std::string Measurement::outfile() const
{
	return _outfile;
}

float Measurement::minCutoff() const
{
	return _minCutoff;
}

float Measurement::maxCutoff() const
{
	return _maxCutoff;
}

bool Measurement::isDeferred() const
{
	return _isDeferred;
}

bool Measurement::isLoadManaged() const
{
	return _isLoadManaged;
}

std::vector<int> Measurement::getLoadManagedStacks() const
{
	return _loadManagedStacks;
}