/**
 * @file Lattice.hpp
 * @author Finn Lasse Buessen
 * @brief Representation of a physical lattice. 
 * 
 * @copyright Copyright (c) 2020
 */

#pragma once
#include <array>
#include <vector>
#include <tuple>
#include "lib/Geometry.hpp"
#include "lib/Assert.hpp"

struct SpinModel;
struct Lattice;
namespace LatticeModelFactory {
	struct LatticeUnitCell;
	struct SpinModelUnitCell;

	std::pair<Lattice *, SpinModel *> newLatticeModel(const LatticeModelFactory::LatticeUnitCell &uc, const LatticeModelFactory::SpinModelUnitCell &spinModelDefinition, const int latticeRange, const std::string &ldfPath);
};

/**
 * @brief Component of a spin operator. 
 */
enum struct SpinComponent : int
{
	X = 0, ///< x-component of a spin. 
	Y = 1, ///< y-component of a spin. 
	Z = 2, ///< z-component of a spin. 
	None = 3 ///< No spin component specified. 
};

/**
 * @brief Structure to describe the portion of the lattice that is relevant in terms of the form sum_j v(i1,j)*v(j,i2). 
 * @details Structure to describe the portion of the lattice that is relevant in terms of the form sum_j v(i1,j)*v(j,i2), 
 * where v(...) is a two-particle vertex function with its two lattice site arguments. 
 * The overlap covers all lattice sites which are within the truncatino range of both vertices. 
 * The LatticeOverlap structure contains a list of representative ids to which the first vertex v(i1,j) is mapped under lattice symmetry reduction for every j within the overlap region. 
 * Similarly, it contains a list of representative ids to which the second vertex v(j,i2) is mapped. 
 * Furthermore, the spin permutations required for the symmetry transformations are stored. 
 * Each tuple of (rid1[i],rid2[i],transformedX1[i],transformedY1[i],transformedZ1[i],transformedX2[i],transformedY2[i],transformedZ2[i]) describes the transformations for a single term in the sum over overlapping lattice sites j. 
 */
struct LatticeOverlap
{
	friend std::pair<Lattice *, SpinModel *> LatticeModelFactory::newLatticeModel(const LatticeModelFactory::LatticeUnitCell &uc, const LatticeModelFactory::SpinModelUnitCell &spinModelDefinition, const int latticeRange, const std::string &ldfPath);
	friend struct Lattice;

public:
	LatticeOverlap() : rid1(nullptr), rid2(nullptr), transformedX1(nullptr), transformedY1(nullptr), transformedZ1(nullptr), transformedX2(nullptr), transformedY2(nullptr), transformedZ2(nullptr), size(0) {}

	/**
	 * @brief Construct a new LatticeOverlap object for a given number of sites. 
	 * 
	 * @param size Number of sites in the overlapping region. 
	 */
	LatticeOverlap(const int size) : size(size)
	{
		rid1 = new int[size];
		rid2 = new int[size];
		transformedX1 = new SpinComponent[size];
		transformedY1 = new SpinComponent[size];
		transformedZ1 = new SpinComponent[size];
		transformedX2 = new SpinComponent[size];
		transformedY2 = new SpinComponent[size];
		transformedZ2 = new SpinComponent[size];
	}

	/**
	 * @brief Construct a new LatticeOverlap object from an existing one.
	 *
	 * @param rhs LatticeOverlap object.
	 */
	LatticeOverlap(const LatticeOverlap &rhs) : LatticeOverlap(rhs.size)
	{
		memcpy(rid1, rhs.rid1, size * sizeof(int));
		memcpy(rid2, rhs.rid2, size * sizeof(int));
		memcpy(transformedX1, rhs.transformedX1, size * sizeof(SpinComponent));
		memcpy(transformedY1, rhs.transformedY1, size * sizeof(SpinComponent));
		memcpy(transformedZ1, rhs.transformedZ1, size * sizeof(SpinComponent));
		memcpy(transformedX2, rhs.transformedX2, size * sizeof(SpinComponent));
		memcpy(transformedY2, rhs.transformedY2, size * sizeof(SpinComponent));
		memcpy(transformedZ2, rhs.transformedZ2, size * sizeof(SpinComponent));
	}

	int *rid1;
	int *rid2;
	SpinComponent *transformedX1;
	SpinComponent *transformedY1;
	SpinComponent *transformedZ1;
	SpinComponent *transformedX2;
	SpinComponent *transformedY2;
	SpinComponent *transformedZ2;
	int size;

protected:
	/**
	 * @brief Destroy the LatticeOverlap object. 
	 */
	~LatticeOverlap()
	{
		delete[] rid1;
		delete[] rid2;
		delete[] transformedX1;
		delete[] transformedY1;
		delete[] transformedZ1;
		delete[] transformedX2;
		delete[] transformedY2;
		delete[] transformedZ2;
	}

	/**
	 * @brief Assignment operator. 
	 * 
	 * @param rhs Object to assign. 
	 * @return LatticeOverlap& Reference to self. 
	 */
	LatticeOverlap &operator=(const LatticeOverlap &rhs)
	{
		size = rhs.size;
		delete[] rid1;
		rid1 = new int[size];
		memcpy(rid1, rhs.rid1, size * sizeof(int));
		delete[] rid2;
		rid2 = new int[size];
		memcpy(rid2, rhs.rid2, size * sizeof(int));
		delete[] transformedX1;
		transformedX1 = new SpinComponent[size];
		memcpy(transformedX1, rhs.transformedX1, size * sizeof(SpinComponent));
		delete[] transformedY1;
		transformedY1 = new SpinComponent[size];
		memcpy(transformedY1, rhs.transformedY1, size * sizeof(SpinComponent));
		delete[] transformedZ1;
		transformedZ1 = new SpinComponent[size];
		memcpy(transformedZ1, rhs.transformedZ1, size * sizeof(SpinComponent));
		delete[] transformedX2;
		transformedX2 = new SpinComponent[size];
		memcpy(transformedX2, rhs.transformedX2, size * sizeof(SpinComponent));
		delete[] transformedY2;
		transformedY2 = new SpinComponent[size];
		memcpy(transformedY2, rhs.transformedY2, size * sizeof(SpinComponent));
		delete[] transformedZ2;
		transformedZ2 = new SpinComponent[size];
		memcpy(transformedZ2, rhs.transformedZ2, size * sizeof(SpinComponent));

		return *this;
	}
};

/**
 * @brief Structure to specify a symmetry-transformed lattice site. 
 * @details The symmetry transformation combines a mapping to a different lattice site and a permutation of spin components.
 * The spin permuation is stored in LatticeSiteDescriptor::spinPermutation, where transformedComponent=spinPermutation[originalComponent].
 */
struct LatticeSiteDescriptor
{
	int rid; ///< Representative id of the transformed lattice site. 
	SpinComponent spinPermutation[3]; ///< Spin permutation involved in the transformation. 
};

/**
 * @brief Lattice iterator object
 */
struct LatticeIterator
{
friend struct Lattice;

public:
	/**
	 * @brief Construct a new LatticeIterator object, which points to the representative id 0 per default. 
	 */
	LatticeIterator() : id(0) {}

	/**
	 * @brief Construct a new LatticeIterator object and initialize it to any representative id. 
	 * 
	 * @param id Representative id to initialize. 
	 */
	LatticeIterator(int id) : id(id) {}

	/**
	 * @brief Comparison operator. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Returns true, if both iterators point to the same representative id. Returns false otherwise. 
	 */
	bool operator==(const LatticeIterator &rhs) const
	{
		return id == rhs.id;
	}

	/**
	 * @brief Negative comparison operator. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return bool Returns false, if both iterators point to the same representative id. Returns true otherwise. 
	 */
	bool operator!=(const LatticeIterator &rhs) const
	{
		return id != rhs.id;
	}

	/**
	 * @brief Subtraction operator. 
	 * 
	 * @param rhs Right hand side operand. 
	 * @return int Difference between iterators. 
	 */
	int operator-(const LatticeIterator &rhs) const
	{
		return id - rhs.id;
	}

	/**
	 * @brief Prefix increment operator. 
	 * 
	 * @return LatticeIterator& Reference to self. 
	 */
	virtual LatticeIterator& operator++()
	{
		++id;
		return *this;
	}

	/**
	 * @brief Prefix decrement operator. 
	 * 
	 * @return LatticeIterator& Reference to self. 
	 */
	virtual LatticeIterator& operator--()
	{
		--id;
		return *this;
	}

protected:
	int id; ///< Representative id which the iterator points to. 
};

/**
 * @brief Sublattice iterator object. 
 */
struct SublatticeIterator : public LatticeIterator
{
public:
	/**
	 * @brief Construct a new SublatticeIterator object over a list of representative ids. 
	 * @details The SublatticeIterator does not take ownership of the list. 
	 * It must be guaranteed that the list remains valid during the entire lifetiem of the iterator object. 
	 * The list mus be terminated by a lattice->end() entry. 
	 * 
	 * @param allowedIds List of representative ids that make up the sublattice. Must be terminated by lattice->end(). 
	 */
	SublatticeIterator(int *allowedIds) : LatticeIterator(*allowedIds), offset(0), allowedIds(allowedIds) {}

	/**
	 * @brief Prefix increment operator. 
	 * 
	 * @return SublatticeIterator& Reference to self. 
	 */
	SublatticeIterator &operator++()
	{
		id = allowedIds[++offset];
		return *this;
	}

	/**
	 * @brief Prefix decrement operator. 
	 * 
	 * @return SublatticeIterator& Reference to self. 
	 */
	SublatticeIterator &operator--()
	{
		id = allowedIds[--offset];
		return *this;
	}

protected:
	int offset; ///< Offset relative to the first entry of SublatticeIterator::allowedIds. 
	int *allowedIds; ///< List of representative ids that make up the sublattice. Memory is not owned by the iterator. 
};

/**
 * @brief Representation of a physical lattice, with symmetry information on two-point correlators. 
 * @details The lattice representation exploits symmetries. A given lattice contains a finite number of sites, which depends on the selected lattice range parameter. 
 * The lattice would store all sites within the specified range (distance measured in units of lattice bonds) around any of the basis sites. 
 *
 * Using lattice symmetries, out of any pair of lattice sites (any two-point correlator), one site can always be mapped to a fixed reference si te at the origin. 
 * The second lattice site, under this mapping, would then be mapped to some other site. 
 * Therefore, the two-point correlators after symmetry reduction only depend on a single lattice site. 
 * The minimum set up these sites, which is required to represent any two-point correlator, is enumerated by the so-called representative id, ranging from zero to some finite value which depends on the lattice geometry. 
 * 
 * The lattice representation provides methods to quickly iterate over the lattice, iterate over certain sublattice regions, or to map any two-point correlator to the corresponding representative site. 
 */
struct Lattice
{
	friend std::pair<Lattice *, SpinModel *> LatticeModelFactory::newLatticeModel(const LatticeModelFactory::LatticeUnitCell &uc, const LatticeModelFactory::SpinModelUnitCell &spinModelDefinition, const int latticeRange, const std::string &ldfPath);

protected:
	/**
	 * @brief Create an uninitialized lattice object. 
	 * @details Constructor should not be called directly. Use LatticeModelFactory::newLatticeModel() to create a new lattice. 
	 */
	Lattice() : size(0), _dataSize(0), _symmetryTable(nullptr), _bufferSites(nullptr), _bufferInvertedSites(nullptr), _bufferOverlapMatrices(nullptr), _bufferBasis(nullptr), _bufferLatticeRange(nullptr) {};

public:
	/**
	 * @brief Destroy the Lattice object. 
	 */
	~Lattice()
	{
		delete[] _bufferBasis;
		for (size_t b = 0; b < _basis.size(); ++b) delete[] _bufferLatticeRange[b];
		delete[] _bufferLatticeRange;
		delete[] _symmetryTable;
		delete[] _bufferSites;
		delete[] _bufferInvertedSites;
		delete[] _bufferOverlapMatrices;
	}

	/**
	 * @brief Retrieve iterator to first lattice site. 
	 * 
	 * @return LatticeIterator Iterator to first lattice site. 
	 */
	LatticeIterator begin() const
	{
		return LatticeIterator(0);
	}

	/**
	 * @brief Retrieve iterator to the lattice site i1=(0,0,0,0).
	 * 
	 * @return LatticeIterator Iterator to the lattice site i1=(0,0,0,0).
	 */
	LatticeIterator zero() const
	{
		return LatticeIterator(0);
	}

	/**
	 * @brief Retrieve iterator to the last+1 lattice site. 
	 * 
	 * @return LatticeIterator Iterator to the last+1 lattice site. 
	 */
	LatticeIterator end() const
	{
		return LatticeIterator(_dataSize);
	}

	/**
	 * @brief Retrieve iterator to the specified representative. 
	 * 
	 * @param rid Representative id of the desired lattice site. 
	 * @return LatticeIterator Iterator to the specified representative site. 
	 */
	LatticeIterator fromParametrization(const int rid) const
	{
		ASSERT(rid < size);
		return LatticeIterator(rid);
	}

	/**
	 * @brief Get the coordinates of a specified lattice site, in units of a1, a2, a3, and b. 
	 * 
	 * @param site Input lattice site. 
	 * @return std::tuple<int, int, int, int> Tuple (n1,n2,n3,nb) which describes the lattice site at coordinate n1*a1+n2*a2+n3*a3+b[nb]. 
	 */
	std::tuple<int, int, int, int> getSiteParameters(const LatticeIterator &site) const
	{
		ASSERT(site - begin() < _dataSize);
		return _geometryTable[site.id];
	}

	/**
	 * @brief Get the position of a lattice site in real space. 
	 * 
	 * @param site Input lattice site. 
	 * @return geometry::Vec3<double> Position of the lattice site in real space. 
	 */
	geometry::Vec3<double> getSitePosition(const LatticeIterator &site) const
	{
		ASSERT(site - begin() < _dataSize);
		return double(std::get<0>(_geometryTable[site.id])) * _bravaisLattice[0] + double(std::get<1>(_geometryTable[site.id])) * _bravaisLattice[1] + double(std::get<2>(_geometryTable[site.id])) * _bravaisLattice[2] + _basis[std::get<3>(_geometryTable[site.id])];
	}

	/**
	 * @brief Transform a pair of lattice sites (i1,i2) and transform it to (0,i2'). 
	 * The transformation may involve a permutation of spin components, whose action will be ignored. 
	 * 
	 * @param i1 Input lattice site. 
	 * @param i2 Input lattice site. 
	 * @return int Representative id of i2'. 
	 */
	int symmetryTransform(const LatticeIterator& i1, const LatticeIterator& i2) const
	{
		ASSERT(_symmetryTable[i1.id * _dataSize + i2.id].rid != -1);
		return _symmetryTable[i1.id * _dataSize + i2.id].rid;
	}

	/**
	 * @brief Transform a pair of lattice sites (i1,i2) and transform it to (0,i2'). 
	 * The transformation may involve a permutation of spin components, in which case the transformation is applied to the spinComponent argument. 
	 * 
	 * @param[in] i1 Input lattice site. 
	 * @param[in] i2 Input lattice site. 
	 * @param[out] spinComponent Spin component to transform under the mapping. 
	 * @return int Representative id of i2'. 
	 */
	int symmetryTransform(const LatticeIterator &i1, const LatticeIterator &i2, SpinComponent &spinComponent) const
	{
		ASSERT(_symmetryTable[i1.id * _dataSize + i2.id].rid != -1);
		if (spinComponent != SpinComponent::None) spinComponent = _symmetryTable[i1.id * _dataSize + i2.id].spinPermutation[static_cast<int>(spinComponent)];
		return _symmetryTable[i1.id * _dataSize + i2.id].rid;
	}

	/**
	 * @brief Transform a pair of lattice sites (i1,i2) and transform it to (0,i2'). 
	 * The transformation may involve a permutation of spin components, in which case the transformation is applied to the spinComponent1 and spinComponent2 arguments. 
	 * 
	 * @param[in] i1 Input lattice site. 
	 * @param[in] i2 Input lattice site. 
	 * @param[out] spinComponent1 Spin component to transform under the mapping. 
	 * @param[out] spinComponent2 Spin component to transform under the mapping. 
	 * @return int Representative id of i2'. 
	 */
	int symmetryTransform(const LatticeIterator &i1, const LatticeIterator &i2, SpinComponent &spinComponent1, SpinComponent &spinComponent2) const
	{
		ASSERT(_symmetryTable[i1.id * _dataSize + i2.id].rid != -1);
		ASSERT(&spinComponent1 != &spinComponent2);

		if (spinComponent1 != SpinComponent::None) spinComponent1 = _symmetryTable[i1.id * _dataSize + i2.id].spinPermutation[static_cast<int>(spinComponent1)];
		if (spinComponent2 != SpinComponent::None) spinComponent2 = _symmetryTable[i1.id * _dataSize + i2.id].spinPermutation[static_cast<int>(spinComponent2)];
		return _symmetryTable[i1.id * _dataSize + i2.id].rid;
	}

	/**
	 * @brief Transform a pair of lattice sites (i1,i2) and transform it to (0,i2'). 
	 * The transformation may involve a permutation of spin components, in which case the transformation is applied to the spinComponent1, spinComponent2, and spinComponent3 arguments. 
	 * 
	 * @param[in] i1 Input lattice site. 
	 * @param[in] i2 Input lattice site. 
	 * @param[out] spinComponent1 Spin component to transform under the mapping. 
	 * @param[out] spinComponent2 Spin component to transform under the mapping. 
	 * @param[out] spinComponent3 Spin component to transform under the mapping. 
	 * @return int Representative id of i2'. 
	 */
	int symmetryTransform(const LatticeIterator &i1, const LatticeIterator &i2, SpinComponent &spinComponent1, SpinComponent &spinComponent2, SpinComponent &spinComponent3) const
	{
		ASSERT(_symmetryTable[i1.id * _dataSize + i2.id].rid != -1);
		ASSERT(&spinComponent1 != &spinComponent2);
		ASSERT(&spinComponent2 != &spinComponent3);
		ASSERT(&spinComponent1 != &spinComponent3);

		if (spinComponent1 != SpinComponent::None) spinComponent1 = _symmetryTable[i1.id * _dataSize + i2.id].spinPermutation[static_cast<int>(spinComponent1)];
		if (spinComponent2 != SpinComponent::None) spinComponent2 = _symmetryTable[i1.id * _dataSize + i2.id].spinPermutation[static_cast<int>(spinComponent2)];
		if (spinComponent3 != SpinComponent::None) spinComponent3 = _symmetryTable[i1.id * _dataSize + i2.id].spinPermutation[static_cast<int>(spinComponent3)];
		return _symmetryTable[i1.id * _dataSize + i2.id].rid;
	}

	/**
	 * @brief Retrieve the lattice overlap of the reference site i1=(0,0,0,0) with some other representative lattice site. 
	 * The overlap is defined as all lattice sites which are within range of both lattice site. 
	 * 
	 * @param rid Representative lattice site to determine overlap with i1. 
	 * @return const LatticeOverlap& Lattice overlap descriptor. 
	 */
	const LatticeOverlap &getOverlap(const int rid) const
	{
		ASSERT(rid < size);
		return _bufferOverlapMatrices[rid];
	}

	/**
	 * @brief List of two-spin correlators (i2,i1), where i1=(0,0,0,0) is the reference site and the list includes all reference sites i2. 
	 * 
	 * @return const LatticeSiteDescriptor* List of symmetry-reduced two-spin correlators (i1,i2). 
	 */
	const LatticeSiteDescriptor *getInvertedSites() const
	{
		return _bufferInvertedSites;
	}

	/**
	 * @brief List of two-spin correlators (i1,i2), where i1=(0,0,0,0) is the reference site and the list includes all reference sites i2. 
	 * 
	 * @return const LatticeSiteDescriptor* List of symmetry-reduced two-spin correlators (i1,i2). 
	 */
	const LatticeSiteDescriptor *getSites() const
	{
		return _bufferSites;
	}

	/**
	 * @brief Retrieve an iterator over all sites which are within range of the lattice site b, where b must describe one of the basis sites. 
	 * @see Lattice::getBasis()
	 * 
	 * @param b Iterator to a basis site.  
	 * @return SublatticeIterator Iterator over all sites in range. 
	 */
	SublatticeIterator getRange(const LatticeIterator &b) const
	{
		ASSERT(std::get<0>(getSiteParameters(b)) == 0);
		ASSERT(std::get<1>(getSiteParameters(b)) == 0);
		ASSERT(std::get<2>(getSiteParameters(b)) == 0);
		ASSERT(std::get<3>(getSiteParameters(b)) < _basis.size());
		
		return SublatticeIterator(_bufferLatticeRange[std::get<3>(getSiteParameters(b))]);
	}

	/**
	 * @brief Retrieve an iterator over all sites which are within range of the lattice site (0,0,0,b). 
	 * 
	 * @param b Basis site index. 
	 * @return SublatticeIterator Iterator over all sites in range. 
	 */
	SublatticeIterator getRange(const int b) const
	{
		ASSERT(b < _basis.size());
		return SublatticeIterator(_bufferLatticeRange[b]);
	}

	/**
	 * @brief Retrieve iterator over all basis sites. 
	 * 
	 * @return SublatticeIterator Iterator over all basis sites. 
	 */
	SublatticeIterator getBasis() const
	{
		return SublatticeIterator(_bufferBasis);
	}

	std::vector<geometry::Vec3<double> > _bravaisLattice; ///< List of the three bravais lattice vectors.  
	std::vector<geometry::Vec3<double> > _basis; ///< List of the positions of all basis sites within a lattice unit cell. 

	int size; ///< Number representative sites. 

protected: 
	std::vector<std::tuple<int, int, int, int>> _geometryTable; ///< Internal storage of real space lattice site positions, stored as tuples (a0, a1, a2, b).

	int _dataSize; ///< Total number of all lattice that we store information about (equivalent to size of Lattice::_geometryTable). 
	
	LatticeSiteDescriptor *_symmetryTable; ///< List of symmetry reductions of spin pairs (id1, id2) linearized as id1*_dataSize+id2. Stores the transformation required to map id1 to zero. 
	LatticeSiteDescriptor *_bufferSites; ///< List of transformed sites (0, rid), for all representatives rid. 
	LatticeSiteDescriptor *_bufferInvertedSites; ///< List of transformed sites (rid, 0), for all representatives rid. 
	LatticeOverlap *_bufferOverlapMatrices; ///< List of lattice overlaps, where the i-th entry is the overlap of the two tuples (0,j)(j,i). 

	int *_bufferBasis; ///< List of representative ids of all basis sites. 
	int **_bufferLatticeRange; ///< Table of lists of all site ids in range of site (0,0,0,b). 
};