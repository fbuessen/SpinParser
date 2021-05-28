/**
 * @file SpinModel.hpp
 * @author Finn Lasse Buessen
 * @brief Representation of a spin model with two-spin interactions. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <tuple>
#include "Lattice.hpp"

/**
 * @brief Spin model representation. 
 * @details A spin model represents the set of two-spin interactions in the Hamiltonian. 
 * The initial conditions of the pf-FRG flow depend on these couplings. 
 * The SpinModel object is therefore passed to the initialization of the EffectiveAction object. 
 */
struct SpinModel
{
public:
	/**
	 * @brief Representation of a two-spin interaction.
	 */
	struct SpinInteraction
	{
		/**
		 * @brief Construct a new SpinInteraction object and initialize all interactions to zero. 
		 */
		SpinInteraction()
		{
			memset(&interactionStrength[0][0], 0, 9 * sizeof(float));
		};

		float interactionStrength[3][3]; ///< Interaction strength, encoded as interactionStrength[s1][s2], where s1 is the x, y, or z (0, 1, or 2) component of the first spin and s2 is the component of the second spin. 
	};

	std::vector<std::pair<LatticeIterator, SpinInteraction> > interactions; ///< List of two-spin interactions in the model. All interactions are between the lattice site pair (s1,s2), where s1 is the reference site and s2 is specified by the LatticeIterator. 
	std::vector<std::string> interactionParameters; ///< List of string-form interaction parameters as used in the task file. The order of the list corresponds to the order in SpinModel::interactions. 
};