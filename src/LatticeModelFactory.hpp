/**
 * @file LatticeModelFactory.hpp
 * @author Finn Lasse Buessen
 * @brief Create lattice representations from a lattice unit cell and specification of spin interactions. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <string>
#include <map>
#include <set>
#include "Lattice.hpp"
#include "SpinModel.hpp"

namespace LatticeModelFactory
{
	/**
	 * @brief Representation of a Lattice site. 
	 */
	struct LatticeSite
	{
	public:
		/**
		 * @brief Construct an uninitialized lattice site object. 
		 */
		LatticeSite();

		/**
		 * @brief Construct a well-defined lattice site. 
		 * 
		 * @param a0 Coordinate in units of the first Bravais lattice vector. 
		 * @param a1 Coordinate in units of the second Bravais lattice vector. 
		 * @param a2 Coordinate in units of the third Bravais lattice vector. 
		 * @param b Basis site index. 
		 */
		LatticeSite(const int a0, const int a1, const int a2, const int b);

		/**
		 * @brief Comparison operator. 
		 * 
		 * @param rhs Right hand side operand. 
		 * @return bool Returns true if both objects desribe the same lattice site. Returns false otherwise. 
		 */
		bool operator==(const LatticeSite& rhs) const;

		/**
		 * @brief Negative comparison operator. 
		 * 
		 * @param rhs Right hand side operand. 
		 * @return bool Returns false if both objects desribe the same lattice site. Returns true otherwise. 
		 */
		bool operator!=(const LatticeSite& rhs) const;

		int a0; ///< Coordinate in units of the first Bravais lattice vector. 
		int a1; ///< Coordinate in units of the second Bravais lattice vector. 
		int a2; ///< Coordinate in units of the third Bravais lattice vector. 
		int b; ///< Basis site index. 
	};

	/**
	 * @brief Representation of a lattice bond. 
	 */
	struct LatticeBond
	{
	public:
		/**
		 * @brief Construct a well-defined lattice bond object. 
		 * 
		 * @param fromB Basis index of the site from which the bond emanates.  
		 * @param toB Basis index of the site which the bond connects to.  
		 * @param da0 Difference of the two connecting sites in units of the first Bravais lattice vector. 
		 * @param da1 Difference of the two connecting sites in units of the second Bravais lattice vector. 
		 * @param da2 Difference of the two connecting sites in units of the third Bravais lattice vector. 
		 */
		LatticeBond(const int fromB, const int toB, const int da0, const int da1, const int da2);

		/**
		 * @brief Determine whether the lattice bond is attached to a specified site. 
		 * 
		 * @param site Site to test. 
		 * @return bool Returns true if one of the ends of the lattice bond is connected the the specified site. Returns false otherwise. 
		 */
		bool isAttachedToSite(const LatticeSite& site) const;

		/**
		 * @brief Determine whether the lattice bond is connecting two specified sites, ignoring the orientation of the bond. 
		 * 
		 * @param site1 First lattice site. 
		 * @param site2 Second lattice site. 
		 * @return bool Returns true if the bond connects the two sites, returns false otherwise. 
		 */
		bool isConnectingSites(const LatticeSite& site1, const LatticeSite& site2) const;

		/**
		 * @brief Determine whether the lattice bond is connecting two specified sites, checking also the orientation of the bond. 
		 * 
		 * @param siteFrom Site from which the bond should emanate. 
		 * @param siteTo Site to which the bond should connect. 
		 * @return bool Returns true if the bonds connects the two sites, returns false otherwise. 
		 */
		bool isConnectingFromTo(const LatticeSite& siteFrom, const LatticeSite& siteTo) const;

		/**
		 * @brief Given one lattice site, determine the other site which the bond is connected to. 
		 * If the site is connecte to the specified lattice sites multiple times (i.e. in a lattice with a monoatomic basis), returns all other connecting sites. 
		 * 
		 * @param site First lattice site. 
		 * @return std::vector<LatticeSite> Other sites which the bond is connected to. 
		 */
		std::vector<LatticeSite> getOtherEnd(const LatticeSite& site) const;

		int fromB; ///< Basis index of the site from which the bond emanates.  
		int toB; ///< Basis index of the site which the bond connects to.  
		int da0; ///<Difference of the two connecting sites in units of the first Bravais lattice vector. 
		int da1; ///<Difference of the two connecting sites in units of the second Bravais lattice vector. 
		int da2; ///<Difference of the two connecting sites in units of the third Bravais lattice vector. 
	};

	/**
	 * @brief Representation of a lattice unit cell. 
	 */
	struct LatticeUnitCell
	{
	public:
		/**
		 * @brief Construct an unitialized lattice unit cell. 
		 */
		LatticeUnitCell();

		/**
		 * @brief Construct a new LatticeUnitCell and initialize it from a specification file placed in the given resource bundle (directory). 
		 * This will search all .xml files in the specified directory for valid unit cell definitions and use the first one that matches the desired name. 
		 * 
		 * @param latticeName Name of the unit cell specification to search for. 
		 * @param bundle Directory to search. 
		 */
		LatticeUnitCell(const std::string &latticeName, const std::string &bundle);

		std::vector<geometry::Vec3<double> > latticeVectors; ///< List of the three Bravais lattice vectors. 
		std::vector<geometry::Vec3<double> > basisSites; ///< List of all basis site positions. 
		std::vector<LatticeBond> latticeBonds; ///< List of all lattice bonds. 

	private:
		/**
		 * @brief Initialize the internal data from a specification file placed in the given resource bundle. 
		 * 
		 * @param latticeName Name of the unit cell specification to search for. 
		 * @param bundle Directory to search. 
		 * @return bool Returns true if a matching lattice unit cell specification was found. Returns false otherwise. 
		 */
		bool _initFromResBundle(const std::string &latticeName, const std::string &bundle);
	};

	/**
	 * @brief Representation of a two-spin interaction. 
	 */
	struct SpinInteraction
	{
	public:
		/**
		 * @brief Construct an uninitialized spin interaction object. 
		 */
		SpinInteraction();

		/**
		 * @brief Construct a two-spin interaction between two lattice sites and set the interaction strength to zero. 
		 *
		 * @param from First interaction site. 
		 * @param to Second interaction site. 
		 */
		SpinInteraction(const LatticeSite& from, const LatticeSite& to);

		/**
		 * @brief Determine whether the spin interaction couples two lattice sites. 
		 * 
		 * @param siteFrom Site from which the interaction emanates. 
		 * @param siteTo Site to which the interaction couples. 
		 * @return bool Returns true if the interaction couples the two lattice sites. Returns false otherwise. 
		 */
		bool isConnectingFromTo(const LatticeSite& siteFrom, const LatticeSite& siteTo) const;

		/**
		 * @brief Determine whether the spin interaction couples two lattice sites, neglecting its orientation. 
		 * 
		 * @param site1 Site from which the interaction emanates.  
		 * @param site2 Site to which the interaction couples. 
		 * @return int Returns 1 if the interaction connects the two sites and the orientation matches. Returns -1 if the orientation does not match. Returns zero if the sites are not connecte by the interaction. 
		 */
		int isConnectingSites(const LatticeSite& site1, const LatticeSite& site2) const;

		/**
		 * @brief Addition assignment operator. 
		 * Adds the interaction strength of a second SpinInteraction object describing coupling between the same two sites. 
		 * 
		 * @param rhs Right hand site operand. 
		 * @return SpinInteraction& Reference to self. 
		 */
		SpinInteraction& operator+=(const SpinInteraction& rhs);

		LatticeSite from; ///< Site from which the interaction emanates. 
		LatticeSite to; ///< Site to which the interaction couples. 
		float interactionStrength[3][3]; ///< Interaction strength, encoded as interactionStrength[s1][s2], where s1 is the x, y, or z (0, 1, or 2) component of the first spin and s2 is the component of the second spin. 
	};

	/**
	 * @brief Spin model unit cell representation. 
	 */
	struct SpinModelUnitCell
	{
	public:
		/**
		 * @brief Construct an uninitialized spin model unit cell. 
		 */
		SpinModelUnitCell();

		/**
		 * @brief Construct a SpinModelUnitCell object and initialize it from a spin model specification with a given name, located in a given resource bundle (directory). 
		 * The coupling strengths are initialized from the modelOptions argument. 
		 * 
		 * @param modelName Name of the model specification to search. 
		 * @param bundle Directory to search for the model specification. 
		 * @param modelOptions List of model options as specified in the task file, to initialize coupling strength. 
		 */
		SpinModelUnitCell(const std::string &modelName, const std::string &bundle, const std::map<std::string, std::string> &modelOptions);

		std::vector<SpinInteraction> interactions; ///< List of spin interactions in the unit cell. 
		std::set<std::string> interactionParameters; ///< List of interaction parameter names as used in the specification file. 

	private:
		/**
		 * @brief Initialize object from a spin model specification with a given name, located in a given resource bundle (directory). 
		 * The coupling strengths are initialized from the modelOptions argument. 
		 * 
		 * @param modelName Name of the model specification to search. 
		 * @param bundle Directory to search for the model specification. 
		 * @param modelOptions List of model options as specified in the task file, to initialize coupling strength. 
		 * @return bool Returns true if a matching unit cell specification was found. Returns false otherwise. 
		 */
		bool _initFromResBundle(const std::string &modelName, const std::string &bundle, const std::map<std::string, std::string> &modelOptions);
	};

	/**
	 * @brief Create and return lattice and spin model objects from given unit cell definitions. 
	 * 
	 * @param uc Lattice unit cell representation. 
	 * @param spinModelDefinition Spin model representation. 
	 * @param latticeRange Lattice range. 
	 * @param ldfPath File name to write the ldf file of the newly created lattice to. If left empty, do not write an ldf file. 
	 * @return std::pair<Lattice *, SpinModel *> Newly created Lattice and SpinModel objects. 
	 */
	std::pair<Lattice *, SpinModel *> newLatticeModel(const LatticeUnitCell &uc, const SpinModelUnitCell &spinModelDefinition, const int latticeRange, const std::string &ldfPath = "");
};