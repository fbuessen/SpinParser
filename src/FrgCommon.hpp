/**
 * @file FrgCommon.hpp
 * @author Finn Lasse Buessen
 * @brief Hub for central objects in pf-FRG calculations. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include "FrequencyDiscretization.hpp"
#include "CutoffDiscretization.hpp"
#include "Lattice.hpp"

/**
 * @brief Hub for central objects in pf-FRG calculations.
 */
struct FrgCommon
{
	friend class SpinParser;
public:
	/**
	 * @brief Retrieve the lattice representation. 
	 * 
	 * @return const Lattice& Lattice representation. 
	 */
	static const Lattice &lattice()
	{
		return *_lattice;
	}

	/**
	 * @brief Retrieve the Matsubara frequency discretization. 
	 * 
	 * @return const FrequencyDiscretization& Matsubara frequency discretization. 
	 */
	static const FrequencyDiscretization &frequency()
	{
		return *_frequency;
	}

	/**
	 * @brief Retrieve the frequency cutoff discretization.
	 * 
	 * @return const CutoffDiscretization& Frequency cutoff discretization. 
	 */
	static const CutoffDiscretization &cutoff()
	{
		return *_cutoff;
	}

private:
	static Lattice *_lattice; ///< Lattice representation. 
	static FrequencyDiscretization *_frequency; ///< Matsubara frequency discretization. 
	static CutoffDiscretization *_cutoff; ///< Frequency cutoff discretization. 
};