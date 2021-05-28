/**
 * @file FrgCommon.cpp
 * @author Finn Lasse Buessen
 * @brief Hub for central objects in pf-FRG calculations. 
 * 
 * @copyright Copyright (c) 2020
 */

#include "FrgCommon.hpp"

Lattice *FrgCommon::_lattice = nullptr;
FrequencyDiscretization *FrgCommon::_frequency = nullptr;
CutoffDiscretization *FrgCommon::_cutoff = nullptr;