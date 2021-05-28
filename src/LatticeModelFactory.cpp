/**
 * @file LatticeModelFactory.cpp
 * @author Finn Lasse Buessen
 * @brief Create lattice representations from a lattice unit cell and specification of spin interactions. 
 * 
 * @copyright Copyright (c) 2020
 */

#include "LatticeModelFactory.hpp"
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/regex.hpp>
#include <boost/format.hpp>
#include "lib/Geometry.hpp"
#include "lib/InputParser.hpp"
#include "lib/Exception.hpp"
#include "lib/Log.hpp"

#define __EPSILON 0.00001
#define PI 3.14159265358979323846

namespace LatticeModelFactory
{
	#pragma region LatticeSite definition
	LatticeSite::LatticeSite()
	{
		a0 = 0;
		a1 = 0;
		a2 = 0;
		b = 0;
	}

	LatticeSite::LatticeSite(const int a0, const int a1, const int a2, const int b) : LatticeSite()
	{
		this->a0 = a0;
		this->a1 = a1;
		this->a2 = a2;
		this->b = b;
	}

	bool LatticeSite::operator==(const LatticeSite& rhs) const
	{
		if (a0 == rhs.a0 && a1 == rhs.a1 && a2 == rhs.a2 && b == rhs.b) return true;
		else return false;
	}

	bool LatticeSite::operator!=(const LatticeSite& rhs) const
	{
		return !operator==(rhs);
	}
	#pragma endregion

	#pragma region LatticeBond definition
	LatticeBond::LatticeBond(const int fromB, const int toB, const int da0, const int da1, const int da2)
	{
		this->fromB = fromB;
		this->toB = toB;
		this->da0 = da0;
		this->da1 = da1;
		this->da2 = da2;
	}

	bool LatticeBond::isAttachedToSite(const LatticeSite& site) const
	{
		if (fromB == site.b || toB == site.b) return true;
		else return false;
	}

	bool LatticeBond::isConnectingSites(const LatticeSite& site1, const LatticeSite& site2) const
	{
		if (isConnectingFromTo(site1, site2)) return true;
		else if (isConnectingFromTo(site2, site1)) return true;
		else return false;
	}

	bool LatticeBond::isConnectingFromTo(const LatticeSite& siteFrom, const LatticeSite& siteTo) const
	{
		if (da0 == (siteTo.a0 - siteFrom.a0) && da1 == (siteTo.a1 - siteFrom.a1) && da2 == (siteTo.a2 - siteFrom.a2) && fromB == siteFrom.b && toB == siteTo.b) return true;
		else return false;
	}

	std::vector<LatticeSite> LatticeBond::getOtherEnd(const LatticeSite& site) const
	{
		if (isAttachedToSite(site))
		{
			if (fromB == site.b)
			{
				if (fromB == toB) return std::vector<LatticeSite>({ LatticeSite(site.a0 + da0, site.a1 + da1, site.a2 + da2, toB), LatticeSite(site.a0 - da0, site.a1 - da1, site.a2 - da2, fromB) });
				else return std::vector<LatticeSite>({ LatticeSite(site.a0 + da0, site.a1 + da1, site.a2 + da2, toB) });
			}
			else return std::vector<LatticeSite>({ LatticeSite(site.a0 - da0, site.a1 - da1, site.a2 - da2, fromB) });
		}
		else throw Exception(Exception::Type::ArgumentError, "Lattice bond is not attached to specified site");
	}
	#pragma endregion

	#pragma region LatticeUniteCell definition
	LatticeUnitCell::LatticeUnitCell() {}

	LatticeUnitCell::LatticeUnitCell(const std::string &latticeName, const std::string &bundle)
	{
		if (!_initFromResBundle(latticeName, bundle)) throw Exception(Exception::Type::InitializationError, "Invalid lattice unit cell definition.");
	}

	bool LatticeUnitCell::_initFromResBundle(const std::string &latticeName, const std::string &bundle)
	{
		boost::filesystem::directory_iterator end;
		for (boost::filesystem::directory_iterator it(bundle); it != end; ++it)
		{
			if (boost::filesystem::is_regular_file(it->path()) && it->path().extension() == ".xml")
			{
				boost::property_tree::ptree latticeDefinition;
				boost::property_tree::xml_parser::read_xml(it->path().string(), latticeDefinition);

				for (auto u : latticeDefinition)
				{
					if (u.first != "unitcell") continue;
					auto uc = u.second;
					if (!uc.get_optional<std::string>("<xmlattr>.name")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell. Unit cell name is undefined.");

					if (uc.get<std::string>("<xmlattr>.name") == latticeName)
					{
						//read unit cell
						if (uc.count("primitive") != 3) throw Exception(Exception::Type::InitializationError, "Invalid unit cell. Unit cell must define three primitive vectors.");
						if (uc.count("site") == 0) throw Exception(Exception::Type::InitializationError, "Invalid unit cell. Unit cell must contain at least one lattice site.");
						if (uc.count("bond") == 0) throw Exception(Exception::Type::InitializationError, "Invalid unit cell. Unit cell must contain at least one lattice bond.");

						for (auto p : uc)
						{
							//lattice vector
							if (p.first == "primitive")
							{
								if (!p.second.get_optional<std::string>("<xmlattr>.x")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell primitive. x attribute missing.");
								if (!p.second.get_optional<std::string>("<xmlattr>.y")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell primitive. y attribute missing.");
								if (!p.second.get_optional<std::string>("<xmlattr>.z")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell primitive. z attribute missing.");

								latticeVectors.push_back(geometry::Vec3<double>(
									InputParser::stringToDouble(p.second.get<std::string>("<xmlattr>.x")),
									InputParser::stringToDouble(p.second.get<std::string>("<xmlattr>.y")),
									InputParser::stringToDouble(p.second.get<std::string>("<xmlattr>.z"))
									));
							}
							//basis site
							else if (p.first == "site")
							{
								if (!p.second.get_optional<std::string>("<xmlattr>.x")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell site. x attribute missing.");
								if (!p.second.get_optional<std::string>("<xmlattr>.y")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell site. y attribute missing.");
								if (!p.second.get_optional<std::string>("<xmlattr>.z")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell site. z attribute missing.");

								basisSites.push_back(geometry::Vec3<double>(
									InputParser::stringToDouble(p.second.get<std::string>("<xmlattr>.x")),
									InputParser::stringToDouble(p.second.get<std::string>("<xmlattr>.y")),
									InputParser::stringToDouble(p.second.get<std::string>("<xmlattr>.z"))
									));
							}
							//bond
							else if (p.first == "bond")
							{
								if (!p.second.get_optional<std::string>("<xmlattr>.from")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell bond. 'from' attribute missing. ");
								if (!p.second.get_optional<std::string>("<xmlattr>.to")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell bond. 'to' attribute missing. ");
								if (!p.second.get_optional<std::string>("<xmlattr>.dx")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell bond. 'dx' attribute missing. ");
								if (!p.second.get_optional<std::string>("<xmlattr>.dy")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell bond. 'dy' attribute missing. ");
								if (!p.second.get_optional<std::string>("<xmlattr>.dz")) throw Exception(Exception::Type::InitializationError, "Invalid unit cell bond. 'dz attribute missing. ");

								int from = std::stoi(p.second.get<std::string>("<xmlattr>.from"));
								int to = std::stoi(p.second.get<std::string>("<xmlattr>.to"));
								int dx = std::stoi(p.second.get<std::string>("<xmlattr>.dx"));
								int dy = std::stoi(p.second.get<std::string>("<xmlattr>.dy"));
								int dz = std::stoi(p.second.get<std::string>("<xmlattr>.dz"));
								this->latticeBonds.push_back(LatticeBond(from, to, dx, dy, dz));
							}
						}
						return true;
					}
				}
			}
		}
		return false;
	}
	#pragma endregion

	#pragma region SpinInteraction definition
	SpinInteraction::SpinInteraction()
	{
		memset(&interactionStrength[0][0], 0, 9 * sizeof(float));
	};

	SpinInteraction::SpinInteraction(const LatticeSite& from, const LatticeSite& to) : from(from), to(to)
	{
		memset(&interactionStrength[0][0], 0, 9 * sizeof(float));
	};

	bool SpinInteraction::isConnectingFromTo(const LatticeSite& siteFrom, const LatticeSite& siteTo) const
	{
		if ((to.a0 - from.a0) == (siteTo.a0 - siteFrom.a0) && (to.a1 - from.a1) == (siteTo.a1 - siteFrom.a1) && (to.a2 - from.a2) == (siteTo.a2 - siteFrom.a2) && from.b == siteFrom.b && to.b == siteTo.b) return true;
		else return false;
	}

	int SpinInteraction::isConnectingSites(const LatticeSite& site1, const LatticeSite& site2) const
	{
		if (isConnectingFromTo(site1, site2)) return 1;
		else if (isConnectingFromTo(site2, site1)) return -1;
		else return 0;
	}

	SpinInteraction& SpinInteraction::operator+=(const SpinInteraction& rhs)
	{
		if (from == rhs.from && to == rhs.to)
		{
			for (int s1 = 0; s1 < 3; ++s1)
			{
				for (int s2 = 0; s2 < 3; ++s2)
				{
					interactionStrength[s1][s2] += rhs.interactionStrength[s1][s2];
				}
			}
		}
		else if (from == rhs.to && to == rhs.from)
		{
			for (int s1 = 0; s1 < 3; ++s1)
			{
				for (int s2 = 0; s2 < 3; ++s2)
				{
					interactionStrength[s1][s2] += rhs.interactionStrength[s2][s1];
				}
			}
		}
		else throw Exception(Exception::Type::ArgumentError, "Cannnot combine incompatible spin interactions.");

		return *this;
	}
	#pragma endregion

	#pragma region SpinModelUnitCell definition
	SpinModelUnitCell::SpinModelUnitCell() {}

	SpinModelUnitCell::SpinModelUnitCell(const std::string &modelName, const std::string &bundle, const std::map<std::string, std::string> &modelOptions)
	{
		if (!_initFromResBundle(modelName, bundle, modelOptions)) throw Exception(Exception::Type::InitializationError, "Invalid spin model definition.");
	}

	bool SpinModelUnitCell::_initFromResBundle(const std::string &res, const std::string &bundle, const std::map<std::string, std::string> &modelOptions)
	{
		//read spin model
		boost::filesystem::directory_iterator end;
		for (boost::filesystem::directory_iterator it(bundle); it != end; ++it)
		{
			if (boost::filesystem::is_regular_file(it->path()) && it->path().extension() == ".xml")
			{
				boost::property_tree::ptree modelDefinition;
				boost::property_tree::xml_parser::read_xml(it->path().string(), modelDefinition);

				for (auto m : modelDefinition)
				{
					if (m.first != "model") continue;
					auto model = m.second;
					if (!model.get_optional<std::string>("<xmlattr>.name")) throw Exception(Exception::Type::InitializationError, "Invalid spin model. Model name is undefined.");

					if (model.get<std::string>("<xmlattr>.name") == res)
					{
						//read model
						if (model.count("interaction") == 0) throw Exception(Exception::Type::InitializationError, "Invalid spin model. Spin model must define at least one interaction.");

						for (auto i : model)
						{
							if (i.first != "interaction") continue;

							if (!i.second.get_optional<std::string>("<xmlattr>.parameter")) throw Exception(Exception::Type::InitializationError, "Invalid interaction. Parameter name not specified. ");
							if (!i.second.get_optional<std::string>("<xmlattr>.type")) throw Exception(Exception::Type::InitializationError, "Invalid interaction. Interaction type not specified. ");
							if (!i.second.get_optional<std::string>("<xmlattr>.from")) throw Exception(Exception::Type::InitializationError, "Invalid interaction. No target sites specified. ");
							if (!i.second.get_optional<std::string>("<xmlattr>.to")) throw Exception(Exception::Type::InitializationError, "Invalid interaction. No target sites specified. ");

							std::string parameter = i.second.get<std::string>("<xmlattr>.parameter");
							std::string from = i.second.get<std::string>("<xmlattr>.from");
							std::string to = i.second.get<std::string>("<xmlattr>.to");
							std::string type = i.second.get<std::string>("<xmlattr>.type");

							//parse lattice sites
							boost::regex pattern("(-{0,1}\\d+)[, ](-{0,1}\\d+)[, ](-{0,1}\\d+)[, ](\\d+)");
							boost::smatch matchSite1;
							boost::smatch matchSite2;
							if (!boost::regex_match(from, matchSite1, pattern)) throw Exception(Exception::Type::InitializationError, "Invalid spin model. Site '" + from + "' is ill-defined. ");
							if (!boost::regex_match(to, matchSite2, pattern)) throw Exception(Exception::Type::InitializationError, "Invalid spin model. Site '" + to + "' is ill-defined. ");

							//parse coupling strength
							float interactionStrength;
							if (modelOptions.count(parameter)) interactionStrength = InputParser::stringToFloat(modelOptions.at(parameter));
							else throw Exception(Exception::Type::InitializationError, "Interaction parameter '" + parameter + "' is not defined in the taskfile");

							//generate interaction
							SpinInteraction interaction(LatticeSite(std::stoi(matchSite1[1].str()), std::stoi(matchSite1[2].str()), std::stoi(matchSite1[3].str()), std::stoi(matchSite1[4].str())), LatticeSite(std::stoi(matchSite2[1].str()), std::stoi(matchSite2[2].str()), std::stoi(matchSite2[3].str()), std::stoi(matchSite2[4].str())));
							if (type == "heisenberg")
							{
								interaction.interactionStrength[0][0] += interactionStrength;
								interaction.interactionStrength[1][1] += interactionStrength;
								interaction.interactionStrength[2][2] += interactionStrength;
							}
							else if (type == "xxyy")
							{
								interaction.interactionStrength[0][0] += interactionStrength;
								interaction.interactionStrength[1][1] += interactionStrength;
							}
							else if (type == "gx")
							{
								interaction.interactionStrength[1][2] += interactionStrength;
								interaction.interactionStrength[2][1] += interactionStrength;
							}
							else if (type == "gy")
							{
								interaction.interactionStrength[2][0] += interactionStrength;
								interaction.interactionStrength[0][2] += interactionStrength;
							}
							else if (type == "gz")
							{
								interaction.interactionStrength[0][1] += interactionStrength;
								interaction.interactionStrength[1][0] += interactionStrength;
							}
							else
							{
								boost::regex pattern("(-?)([xyz])([xyz])");
								boost::smatch match;
								if (boost::regex_match(type, match, pattern))
								{
									float sign = (match[1] == "-") ? -1.0f : 1.0f;
									int s1 = 0;
									int s2 = 0;
									if (match[2] == "x") s1 = 0;
									else if (match[2] == "y") s1 = 1;
									else if (match[2] == "z") s1 = 2;
									if (match[3] == "x") s2 = 0;
									else if (match[3] == "y") s2 = 1;
									else if (match[3] == "z") s2 = 2;

									interaction.interactionStrength[s1][s2] += sign * interactionStrength;
								}
								else throw Exception(Exception::Type::InitializationError, "Invalid spin model. Unknown interaction type '" + type + "'.");
							}

							//add interaction to spin model
							auto addInteraction = [&](const SpinInteraction &i)->void
							{
								for (auto j = interactions.begin(); j != interactions.end(); ++j)
								{
									if (j->isConnectingSites(i.from, i.to) != 0)
									{
										(*j) += i;
										return;
									}
								}
								interactions.push_back(i);
							};
							addInteraction(interaction);

							//add required interaction parameter
							interactionParameters.insert(parameter);
						}
						return true;
					}
				}
			}
		}
		return false;
	};
	#pragma endregion

	#pragma region private interface
	struct SpinPermutation
	{
	public:
		SpinPermutation()
		{
			transformedComponent[0] = SpinComponent::X;
			transformedComponent[1] = SpinComponent::Y;
			transformedComponent[2] = SpinComponent::Z;
		}

		SpinPermutation(SpinComponent transformedX, SpinComponent transformedY, SpinComponent transformedZ)
		{
			transformedComponent[0] = transformedX;
			transformedComponent[1] = transformedY;
			transformedComponent[2] = transformedZ;
		}

		SpinPermutation(const int n)
		{
			switch (n)
			{
			case 0:
				transformedComponent[0] = SpinComponent::X;
				transformedComponent[1] = SpinComponent::Y;
				transformedComponent[2] = SpinComponent::Z;
				break;
			case 1:
				transformedComponent[0] = SpinComponent::X;
				transformedComponent[1] = SpinComponent::Z;
				transformedComponent[2] = SpinComponent::Y;
				break;
			case 2:
				transformedComponent[0] = SpinComponent::Y;
				transformedComponent[1] = SpinComponent::X;
				transformedComponent[2] = SpinComponent::Z;
				break;
			case 3:
				transformedComponent[0] = SpinComponent::Y;
				transformedComponent[1] = SpinComponent::Z;
				transformedComponent[2] = SpinComponent::X;
				break;
			case 4:
				transformedComponent[0] = SpinComponent::Z;
				transformedComponent[1] = SpinComponent::X;
				transformedComponent[2] = SpinComponent::Y;
				break;
			case 5:
				transformedComponent[0] = SpinComponent::Z;
				transformedComponent[1] = SpinComponent::Y;
				transformedComponent[2] = SpinComponent::X;
				break;
			default:
				throw Exception(Exception::Type::ArgumentError, "Specified spin permutation does not exist");
				break;
			}
		}

		static SpinPermutation identity()
		{
			return SpinPermutation(SpinComponent::X, SpinComponent::Y, SpinComponent::Z);
		}

		static SpinPermutation inverse(const SpinPermutation& rhs)
		{
			SpinPermutation p;
			if (rhs.transformedComponent[0] == SpinComponent::X) p.transformedComponent[0] = SpinComponent::X;
			else if (rhs.transformedComponent[0] == SpinComponent::Y) p.transformedComponent[1] = SpinComponent::X;
			else if (rhs.transformedComponent[0] == SpinComponent::Z) p.transformedComponent[2] = SpinComponent::X;

			if (rhs.transformedComponent[1] == SpinComponent::X) p.transformedComponent[0] = SpinComponent::Y;
			else if (rhs.transformedComponent[1] == SpinComponent::Y) p.transformedComponent[1] = SpinComponent::Y;
			else if (rhs.transformedComponent[1] == SpinComponent::Z) p.transformedComponent[2] = SpinComponent::Y;

			if (rhs.transformedComponent[2] == SpinComponent::X) p.transformedComponent[0] = SpinComponent::Z;
			else if (rhs.transformedComponent[2] == SpinComponent::Y) p.transformedComponent[1] = SpinComponent::Z;
			else if (rhs.transformedComponent[2] == SpinComponent::Z) p.transformedComponent[2] = SpinComponent::Z;

			return p;
		}

		bool operator==(const SpinPermutation& rhs) const
		{
			return transformedComponent[0] == rhs.transformedComponent[0] && transformedComponent[1] == rhs.transformedComponent[1] && transformedComponent[2] == rhs.transformedComponent[2];
		}

		SpinComponent transformedComponent[3];
	};

	//return a list of all nearest neighbors of given lattice site
	std::vector<LatticeSite> getNeighbors(const LatticeUnitCell& uc, const LatticeSite& site)
	{
		std::vector<LatticeSite> neighbors;

		for (auto bond : uc.latticeBonds)
		{
			if (bond.isAttachedToSite(site))
			{
				std::vector<LatticeSite> newNeighbors = bond.getOtherEnd(site);
				neighbors.insert(neighbors.end(), newNeighbors.begin(), newNeighbors.end());
			}
		}
		return neighbors;
	}

	//return a list of all neighbors of given lattice site within a specified range
	std::vector<LatticeSite> constructRangeAroundSite(const LatticeUnitCell& uc, const LatticeSite& site, int range)
	{
		std::vector<LatticeSite> sites({ site });
		for (int i = 0; i < range; ++i)
		{
			//for each site of the set, add all neighbors
			std::vector<LatticeSite> neighbors;
			for (auto site : sites)
			{
				std::vector<LatticeSite> ns = getNeighbors(uc, site);
				for (auto n : ns) if (std::find(neighbors.begin(), neighbors.end(), n) == neighbors.end()) neighbors.push_back(n);
			}
			for (auto n : neighbors) if (std::find(sites.begin(), sites.end(), n) == sites.end()) sites.push_back(n);
		}
		return sites;
	}

	//return the coordinates of given lattice site
	geometry::Vec3<double> getSitePosition(const LatticeUnitCell& uc, const LatticeSite& site)
	{
		return double(site.a0) * uc.latticeVectors[0] + double(site.a1) * uc.latticeVectors[1] + double(site.a2) * uc.latticeVectors[2] + uc.basisSites[site.b];
	}

	//return true and write result to LatticeSite &site if a lattice site exists at the specified coordinates, otherwise return false
	bool siteAtPosition(const LatticeUnitCell& uc, const geometry::Vec3<double>& position, LatticeSite& site)
	{
		for (int b = 0; b < int(uc.basisSites.size()); ++b)
		{
			geometry::Vec3<double> n = geometry::Mat3<double>(uc.latticeVectors[0], uc.latticeVectors[1], uc.latticeVectors[2]).inverse() * (position - uc.basisSites[b]);
			if (fabs(n.x - std::round(n.x)) < __EPSILON && fabs(n.y - std::round(n.y)) < __EPSILON && fabs(n.z - std::round(n.z)) < __EPSILON)
			{
				site.a0 = std::lround(n.x);
				site.a1 = std::lround(n.y);
				site.a2 = std::lround(n.z);
				site.b = b;
				return true;
			}
		}
		return false;
	};

	//searches lattice automorphisms f_i that map site2 onto (0,0,0,0). Return f_i(site1). If reducedSearch==true, return after the first f_i has been found. 
	std::vector<std::pair<LatticeSite, SpinPermutation>> symmetryReduce(const LatticeUnitCell& uc, const SpinModelUnitCell& spinModel, const LatticeSite& site1, const LatticeSite& site2, bool reducedSearch = false)
	{
		//define reference site ref
		LatticeSite ref = LatticeSite(0, 0, 0, 0);

		//collect minimal set of lattice sites which needs to be considered in order to verify symmetry transformations
		std::vector<LatticeSite> preImageSites;
		for (unsigned int b = 0; b < uc.basisSites.size(); ++b)
		{
			LatticeSite site(0, 0, 0, b);
			if (std::find(preImageSites.begin(), preImageSites.end(), site) == preImageSites.end()) preImageSites.push_back(site);
		}
		for (auto bond : uc.latticeBonds)
		{
			LatticeSite site(bond.da0, bond.da1, bond.da2, bond.toB);
			if (std::find(preImageSites.begin(), preImageSites.end(), site) == preImageSites.end()) preImageSites.push_back(site);
		}

		//collect minimal set of spin interactions which need to be considered in order to verify symmetry transformations
		std::vector<SpinInteraction> preImageInteractions = spinModel.interactions; //TODO select only those that are within neighborhood[site2.b]

		//choose two neighbors (n1, n2) of site2 such that the sites are non-collinear
		LatticeSite n1, n2;
		std::vector<LatticeSite> neighbors = getNeighbors(uc, site2);
		if (neighbors.size() < 3) throw Exception(Exception::Type::InitializationError, "Could not build lattice. Too few neighbors. ");
		n1 = neighbors[0];
		for (unsigned int i = 1; i < neighbors.size(); ++i)
		{
			if (geometry::cross(getSitePosition(uc, n1) - getSitePosition(uc, site2), getSitePosition(uc, neighbors[i]) - getSitePosition(uc, site2)).norm() > __EPSILON)
			{
				n2 = neighbors[i];
				break;
			}
			if (i == neighbors.size() - 1) throw Exception(Exception::Type::InitializationError, "Could not build lattice. Unfavorable lattice geometry. ");
		}

		//container to store all valid transformations f_i
		std::vector<std::pair<geometry::Mat4<double>, SpinPermutation>> validTransformations;

		//try to match (site2, n1, n2) with any combination of neighbors (ref, c1, c2) where c1 and c2 are neighbors of ref
		std::vector<LatticeSite> refNeighbors = getNeighbors(uc, ref);
		for (auto c1 : refNeighbors)
		{
			for (auto c2 : refNeighbors)
			{
				//skip if (c1, c2, ref) collinear
				if (cross(getSitePosition(uc, c1) - getSitePosition(uc, ref), getSitePosition(uc, c2) - getSitePosition(uc, ref)).norm() < __EPSILON) continue;

				for (int inversion = 0; inversion <= 1; ++inversion)
				{
					if (reducedSearch && validTransformations.size() > 0) break;

					//init transformation
					geometry::Mat4<double> transformation = (inversion) ? geometry::Mat4<double>::inversion() : geometry::Mat4<double>::identity();
					SpinPermutation spinTransformation;

					//move site2 to ref
					transformation = geometry::Mat4<double>::translation(getSitePosition(uc, ref) - transformation * getSitePosition(uc, site2)) * transformation;

					//rotate n1 onto c1
					geometry::Vec3<double> a = (transformation * getSitePosition(uc, n1) - getSitePosition(uc, ref)).normalize();
					geometry::Vec3<double> b = (getSitePosition(uc, c1) - getSitePosition(uc, ref)).normalize();
					geometry::Vec3<double> axis = cross(a, b);
					double angle = acos(dot(a, b));
					if (axis.norm() < __EPSILON)
					{
						//a and b are already antiparallel or parallel. If they are antiparallel, rotate them. Otherwise do nothing. 
						if (dot(a, b) < 0)
						{
							axis = (a.x != 0) ? geometry::Vec3<double>(-a.y, a.x, 0) : geometry::Vec3<double>(0, -a.z, a.y);
							angle = PI;
						}
						else
						{
							axis = geometry::Vec3<double>(1, 0, 0);
							angle = 0;
						}
					}
					transformation = geometry::Mat4<double>::rotation(axis, getSitePosition(uc, ref), angle) * transformation;

					//try to rotate n2 onto c2
					geometry::Vec3<double> x = (transformation * getSitePosition(uc, n2) - getSitePosition(uc, ref) - dot(transformation * getSitePosition(uc, n2) - getSitePosition(uc, ref), (getSitePosition(uc, c1) - getSitePosition(uc, ref)).normalize()) * (getSitePosition(uc, c1) - getSitePosition(uc, ref)).normalize()).normalize();
					geometry::Vec3<double> y = (getSitePosition(uc, c2) - getSitePosition(uc, ref) - dot(getSitePosition(uc, c2) - getSitePosition(uc, ref), (getSitePosition(uc, c1) - getSitePosition(uc, ref)).normalize()) * (getSitePosition(uc, c1) - getSitePosition(uc, ref)).normalize()).normalize();
					axis = cross(x, y);
					angle = acos(dot(x, y));
					if (axis.norm() < __EPSILON)
					{
						axis = b;
						angle = (dot(x, y) < 0) ? PI : 0.0;
					}
					transformation = geometry::Mat4<double>::rotation(axis, getSitePosition(uc, ref), angle) * transformation;

					//perform sanity checks on the basic transformation properties
					//T(site2) should match ref
					double e = (transformation * getSitePosition(uc, site2) - getSitePosition(uc, ref)).norm();
					if (e > __EPSILON) throw Exception(Exception::Type::InternalError, "Lattice symmetry calculation has failed. (Deviation of symmetry transformed site [s2] from target [ref] is " + std::to_string(e) + ", should be zero)");

					//T(n1) should match of c1
					e = (transformation * getSitePosition(uc, n1) - getSitePosition(uc, c1)).norm();
					if (e > __EPSILON) throw Exception(Exception::Type::InternalError, "Lattice symmetry calculation has failed. (Deviation of symmetry transformed site [n1] from target [c1] is " + std::to_string(e) + ", should be zero)");

					//T(n2) should be coplanar with c1-ref and c2-ref
					e = dot(transformation * getSitePosition(uc, n2) - getSitePosition(uc, ref), cross(getSitePosition(uc, c1) - getSitePosition(uc, ref), getSitePosition(uc, c2) - getSitePosition(uc, ref)));
					if (fabs(e) > __EPSILON) throw Exception(Exception::Type::InternalError, "Lattice symmetry calculation has failed. (Symmetry transformed site [n2] is not coplanar with [c1-ref] and [c2-ref]. Deviation is " + std::to_string(e) + ", should be zero)");

					//check if the transformation is valid
					bool validTransformation = true;

					//check if lattice sites are invariant under transformation
					LatticeSite imageSite;
					for (auto p : preImageSites)
					{
						if (!siteAtPosition(uc, transformation * getSitePosition(uc, p), imageSite))
						{
							validTransformation = false;
							break;
						}
					}
					if (!validTransformation) continue;

					//check if spin couplings are invariant under transformation up to a global permutation of spin components
					bool spinTransformationExists = false;
					for (int i = 0; i < 6; ++i)
					{
						SpinPermutation permutation(i);
						bool validPermutation = true;

						//check whether permutation is compatible with all interactions in the model
						for (auto interaction : preImageInteractions)
						{
							//determine transformed interaction
							LatticeSite fromTransformed;
							siteAtPosition(uc, transformation * getSitePosition(uc, interaction.from), fromTransformed);
							LatticeSite toTransformed;
							siteAtPosition(uc, transformation * getSitePosition(uc, interaction.to), toTransformed);

							float targetInteractionStrength[3][3];
							for (int s1 = 0; s1 < 3; ++s1)
							{
								for (int s2 = 0; s2 < 3; ++s2)
								{
									targetInteractionStrength[s1][s2] = 0.0f;
								}
							}

							for (auto i : spinModel.interactions)
							{
								int connection = i.isConnectingSites(fromTransformed, toTransformed);

								if (connection == 1)
								{
									for (int s1 = 0; s1 < 3; ++s1)
									{
										for (int s2 = 0; s2 < 3; ++s2)
										{
											targetInteractionStrength[s1][s2] = i.interactionStrength[s1][s2];
										}
									}
								}
								else if (connection == -1)
								{
									for (int s1 = 0; s1 < 3; ++s1)
									{
										for (int s2 = 0; s2 < 3; ++s2)
										{
											targetInteractionStrength[s1][s2] = i.interactionStrength[s2][s1];
										}
									}
								}
							}

							//verify transformation
							for (int s1 = 0; s1 < 3; ++s1)
							{
								for (int s2 = 0; s2 < 3; ++s2)
								{
									if (interaction.interactionStrength[s1][s2] != targetInteractionStrength[static_cast<int>(permutation.transformedComponent[s1])][static_cast<int>(permutation.transformedComponent[s2])]) validPermutation = false;
								}
							}
						}

						if (validPermutation)
						{
							spinTransformationExists = true;
							spinTransformation = permutation;
							break;
						}
					}
					if (!spinTransformationExists) validTransformation = false;

					if (!validTransformation) continue;
					else validTransformations.push_back(std::pair<geometry::Mat4<double>, SpinPermutation>(transformation, spinTransformation));
				}
			}
		}

		//compute f_i(site1) and return
		std::vector<std::pair<LatticeSite, SpinPermutation>> fsite1;
		if (validTransformations.size() > 0)
		{
			fsite1.reserve(validTransformations.size());

			//build equivalence class from all transformed sites
			for (auto t : validTransformations)
			{
				LatticeSite ts1;
				if (!siteAtPosition(uc, t.first * getSitePosition(uc, site1), ts1)) throw Exception(Exception::Type::InternalError, "Lattice symmetry calculation has failed. (Internal error. Unexpected outcome of an allegedly valid symmetry)");
				std::pair<LatticeSite, SpinPermutation> transformation(ts1, t.second);

				if (std::find(fsite1.begin(), fsite1.end(), transformation) == fsite1.end()) fsite1.push_back(transformation);
			}
		}
		else throw Exception(Exception::Type::InternalError, "Lattice symmetry calculation has failed. (Internal error. Could not establish identity as a valid symmetry operation)");
		return fsite1;
	}
	#pragma endregion

	std::pair<Lattice *, SpinModel *> newLatticeModel(const LatticeUnitCell &uc, const SpinModelUnitCell &spinModelDefinition, const int latticeRange, const std::string &ldfPath)
	{
		Log::log << Log::LogLevel::Info << "Building lattice spin model..." << Log::endl;

		//generate empty lattice 
		Lattice *lattice = new Lattice;

		//construct neighborhoods around basis sites
		std::vector<std::vector<LatticeSite> > neighborhoods;
		for (int b = 0; b < int(uc.basisSites.size()); ++b) neighborhoods.push_back(constructRangeAroundSite(uc, LatticeSite(0, 0, 0, b), latticeRange));

		//set lattice->_bravaisLattice
		lattice->_bravaisLattice = uc.latticeVectors;

		//set lattice->_basis
		lattice->_basis = uc.basisSites;

		//construct lattice parametrization
		Log::log << Log::LogLevel::Info << "\t...finding lattice parametrization" << Log::endl;
		std::vector<std::pair<int, SpinPermutation>> equivalenceClasses;
		equivalenceClasses.resize(neighborhoods[0].size());
		for (int i = 0; i < int(equivalenceClasses.size()); ++i) equivalenceClasses[i] = std::pair<int, SpinPermutation>(-1, SpinPermutation());
		int rid = 0;
		for (int i = 0; i < int(equivalenceClasses.size()); ++i)
		{
			//skip if the equivalence class of this site has already been defined
			if (equivalenceClasses[i].first != -1) continue;

			//compute symmetry related sites
			auto equiv = symmetryReduce(uc, spinModelDefinition, neighborhoods[0][i], LatticeSite(0, 0, 0, 0));

			//define trivial representative of the group
			equivalenceClasses[i].first = rid;
			equivalenceClasses[i].second = SpinPermutation::identity();

			//assign the same representative id (rid) to all equivalent sites
			for (auto e : equiv)
			{
				//determine symmetry related partner, skip if it has been addressed before
				int j = int(std::find(neighborhoods[0].begin(), neighborhoods[0].end(), e.first) - neighborhoods[0].begin());
				if (equivalenceClasses[j].first != -1) continue;

				//e contains a symmetry transformation of site(i) to site(j), i.e. we need to store the inverse transformation for site(j)
				equivalenceClasses[j].first = rid;
				equivalenceClasses[j].second = SpinPermutation::inverse(e.second);
			}
			++rid;
		}
		//sort lattice parametrization to store the trivial representative of each equivalence class upfront
		for (int r = 0; r < rid; ++r)
		{
			for (int i = 0; i < int(equivalenceClasses.size()); ++i)
			{
				if (equivalenceClasses[i].first == r && equivalenceClasses[i].second == SpinPermutation::identity())
				{
					std::swap(equivalenceClasses[r], equivalenceClasses[i]);
					std::swap(neighborhoods[0][r], neighborhoods[0][i]);
					continue;
				}
			}
		}

		//set lattice->size
		lattice->size = rid;

		//generate lattice sites
		std::vector<LatticeSite> sites = neighborhoods[0];
		for (int b = 0; b < int(uc.basisSites.size()); ++b)
		{
			for (auto n : neighborhoods[b]) if (std::find(sites.begin(), sites.end(), n) == sites.end()) sites.push_back(n);
		}

		//set lattice->_geometryTable
		Log::log << Log::LogLevel::Info << "\t...initializing lattice geometry buffers" << Log::endl;
		lattice->_geometryTable.resize(sites.size());
		for (int i = 0; i < int(sites.size()); ++i) lattice->_geometryTable[i] = std::tuple<int, int, int, int>(sites[i].a0, sites[i].a1, sites[i].a2, sites[i].b);

		//set lattice->_dataSize
		lattice->_dataSize = int(sites.size());

		//set lattice->_bufferBasis
		lattice->_bufferBasis = new int[uc.basisSites.size() + 1];
		for (int b = 0; b < int(uc.basisSites.size()); ++b)
		{
			for (int i = 0; i < int(sites.size()); ++i) if (sites[i] == LatticeSite(0, 0, 0, b)) lattice->_bufferBasis[b] = i;
		}
		lattice->_bufferBasis[uc.basisSites.size()] = int(sites.size());

		//set lattice->_bufferLatticeRange
		lattice->_bufferLatticeRange = new int*[uc.basisSites.size()];
		for (int b = 0; b < int(uc.basisSites.size()); ++b)
		{
			lattice->_bufferLatticeRange[b] = new int[neighborhoods[b].size() + 1];
			for (int n = 0; n < int(neighborhoods[b].size()); ++n)
			{
				for (int i = 0; i < int(sites.size()); ++i) if (neighborhoods[b][n] == sites[i]) lattice->_bufferLatticeRange[b][n] = i;
			}
			lattice->_bufferLatticeRange[b][neighborhoods[b].size()] = int(sites.size());
		}

		//generate lattice->_symmetryTable
		Log::log << Log::LogLevel::Info << "\t...calculating lattice symmetries" << Log::endl;
		lattice->_symmetryTable = new LatticeSiteDescriptor[sites.size() * sites.size()];
		for (int i = 0; i < int(sites.size() * sites.size()); ++i)
		{
			lattice->_symmetryTable[i].rid = -1;
			lattice->_symmetryTable[i].spinPermutation[0] = SpinComponent::X;
			lattice->_symmetryTable[i].spinPermutation[1] = SpinComponent::Y;
			lattice->_symmetryTable[i].spinPermutation[2] = SpinComponent::Z;
		}
		//init entries related to lattice parametrization (transformations with site1=LatticeSite(0,0,0,0))
		for (int n = 0; n < int(neighborhoods[0].size()); ++n)
		{
			lattice->_symmetryTable[0 * sites.size() + n].rid = equivalenceClasses[n].first;
			lattice->_symmetryTable[0 * sites.size() + n].spinPermutation[0] = equivalenceClasses[n].second.transformedComponent[0];
			lattice->_symmetryTable[0 * sites.size() + n].spinPermutation[1] = equivalenceClasses[n].second.transformedComponent[1];
			lattice->_symmetryTable[0 * sites.size() + n].spinPermutation[2] = equivalenceClasses[n].second.transformedComponent[2];
		}
		//init entries related to fundamental transformations with site1=LatticeSite(0,0,0,b) for b non-zero)
		for (int b = 1; b < int(uc.basisSites.size()); ++b)
		{
			for (auto n : neighborhoods[b])
			{
				auto equiv = symmetryReduce(uc, spinModelDefinition, n, LatticeSite(0, 0, 0, b), true);
				if (equiv.size() == 0) throw Exception(Exception::Type::InitializationError, "Could not build lattice. Could not find enough symmetries");
				int bid = int(std::find(sites.begin(), sites.end(), LatticeSite(0, 0, 0, b)) - sites.begin());
				int nid = int(std::find(sites.begin(), sites.end(), n) - sites.begin());
				int eid = int(std::find(sites.begin(), sites.end(), equiv[0].first) - sites.begin());

				lattice->_symmetryTable[bid * sites.size() + nid].rid = lattice->_symmetryTable[0 * sites.size() + eid].rid;
				lattice->_symmetryTable[bid * sites.size() + nid].spinPermutation[0] = equivalenceClasses[eid].second.transformedComponent[static_cast<int>(equiv[0].second.transformedComponent[static_cast<int>(SpinComponent::X)])];
				lattice->_symmetryTable[bid * sites.size() + nid].spinPermutation[1] = equivalenceClasses[eid].second.transformedComponent[static_cast<int>(equiv[0].second.transformedComponent[static_cast<int>(SpinComponent::Y)])];
				lattice->_symmetryTable[bid * sites.size() + nid].spinPermutation[2] = equivalenceClasses[eid].second.transformedComponent[static_cast<int>(equiv[0].second.transformedComponent[static_cast<int>(SpinComponent::Z)])];
			}
		}
		//init remaining entries of all overlapping pairs of sites
		for (int s1 = 0; s1 < int(sites.size()); ++s1)
		{
			for (int s2 = 0; s2 < int(sites.size()); ++s2)
			{
				LatticeSite s1p(0, 0, 0, sites[s1].b);
				LatticeSite s2p(sites[s2].a0 - sites[s1].a0, sites[s2].a1 - sites[s1].a1, sites[s2].a2 - sites[s1].a2, sites[s2].b);
				if (std::find(neighborhoods[s1p.b].begin(), neighborhoods[s1p.b].end(), s2p) == neighborhoods[s1p.b].end()) continue;

				int s1pid = int(std::find(sites.begin(), sites.end(), s1p) - sites.begin());
				int s2pid = int(std::find(sites.begin(), sites.end(), s2p) - sites.begin());
				lattice->_symmetryTable[s1 * sites.size() + s2] = lattice->_symmetryTable[s1pid * sites.size() + s2pid];
			}
		}

		//generate lattice->_bufferSites
		lattice->_bufferSites = new LatticeSiteDescriptor[lattice->size];
		for (int i = 0; i < lattice->size; ++i)
		{
			lattice->_bufferSites[i] = lattice->_symmetryTable[0 * sites.size() + i];
		}

		//generate lattice->_bufferInvertedSites
		lattice->_bufferInvertedSites = new LatticeSiteDescriptor[lattice->size];
		for (int i = 0; i < lattice->size; ++i)
		{
			lattice->_bufferInvertedSites[i] = lattice->_symmetryTable[i * sites.size() + 0];
		}

		//generate lattice->_bufferOverlapMatrices
		auto bufferNewOverlapTable = [&](int rid)
		{
			//overlap table buffer
			std::vector<int> overlapRid1;
			std::vector<int> overlapRid2;
			std::vector<SpinComponent> overlapTX1;
			std::vector<SpinComponent> overlapTY1;
			std::vector<SpinComponent> overlapTZ1;
			std::vector<SpinComponent> overlapTX2;
			std::vector<SpinComponent> overlapTY2;
			std::vector<SpinComponent> overlapTZ2;

			//generate overlap buffer
			LatticeSite i1 = sites[0];
			LatticeSite i2 = sites[rid];

			for (int j = 0; j < int(sites.size()); ++j)
			{
				//verify that j is in range of i1
				if (std::find(neighborhoods[i1.b].begin(), neighborhoods[i1.b].end(), sites[j]) == neighborhoods[i1.b].end()) continue;

				//verify that j is in range of i2
				if (std::find(neighborhoods[i2.b].begin(), neighborhoods[i2.b].end(), LatticeSite(sites[j].a0 - i2.a0, sites[j].a1 - i2.a1, sites[j].a2 - i2.a2, sites[j].b)) == neighborhoods[i2.b].end()) continue;

				//symmetry transform the pair (0,j) to obtain rid1 and (j,rid) to obtain rid2
				LatticeSiteDescriptor t1 = lattice->_symmetryTable[0 * sites.size() + j];
				LatticeSiteDescriptor t2 = lattice->_symmetryTable[j * sites.size() + rid];

				overlapRid1.push_back(t1.rid);
				overlapTX1.push_back(t1.spinPermutation[0]);
				overlapTY1.push_back(t1.spinPermutation[1]);
				overlapTZ1.push_back(t1.spinPermutation[2]);
				overlapRid2.push_back(t2.rid);
				overlapTX2.push_back(t2.spinPermutation[0]);
				overlapTY2.push_back(t2.spinPermutation[1]);
				overlapTZ2.push_back(t2.spinPermutation[2]);
			}

			LatticeOverlap overlap(int(overlapRid1.size()));
			memcpy(overlap.rid1, overlapRid1.data(), overlapRid1.size() * sizeof(int));
			memcpy(overlap.rid2, overlapRid2.data(), overlapRid2.size() * sizeof(int));
			memcpy(overlap.transformedX1, overlapTX1.data(), overlapTX1.size() * sizeof(SpinComponent));
			memcpy(overlap.transformedY1, overlapTY1.data(), overlapTY1.size() * sizeof(SpinComponent));
			memcpy(overlap.transformedZ1, overlapTZ1.data(), overlapTZ1.size() * sizeof(SpinComponent));
			memcpy(overlap.transformedX2, overlapTX2.data(), overlapTX2.size() * sizeof(SpinComponent));
			memcpy(overlap.transformedY2, overlapTY2.data(), overlapTY2.size() * sizeof(SpinComponent));
			memcpy(overlap.transformedZ2, overlapTZ2.data(), overlapTZ2.size() * sizeof(SpinComponent));

			return overlap;
		};
		lattice->_bufferOverlapMatrices = new LatticeOverlap[lattice->size];
		for (int rid = 0; rid < lattice->size; ++rid) lattice->_bufferOverlapMatrices[rid] = bufferNewOverlapTable(rid);

		//init SpinModel
		SpinModel* spinModel = new SpinModel();

		//add interaction parameters
		spinModel->interactionParameters.insert(spinModel->interactionParameters.begin(), spinModelDefinition.interactionParameters.begin(), spinModelDefinition.interactionParameters.end());

		//add interaction strengths
		for (int rid = 0; rid < lattice->size; ++rid)
		{
			auto i1 = lattice->getSiteParameters(lattice->fromParametrization(rid));

			for (auto interaction : spinModelDefinition.interactions)
			{
				int connection = interaction.isConnectingSites(LatticeSite(0, 0, 0, 0), LatticeSite(std::get<0>(i1), std::get<1>(i1), std::get<2>(i1), std::get<3>(i1)));

				if (connection == 1)
				{
					SpinModel::SpinInteraction i;
					for (int s1 = 0; s1 < 3; ++s1)
					{
						for (int s2 = 0; s2 < 3; ++s2) i.interactionStrength[s1][s2] = interaction.interactionStrength[s1][s2];
					}
					spinModel->interactions.push_back(std::pair<LatticeIterator, SpinModel::SpinInteraction>(lattice->fromParametrization(rid), i));
				}
				else if (connection == -1)
				{
					SpinModel::SpinInteraction i;
					for (int s1 = 0; s1 < 3; ++s1)
					{
						for (int s2 = 0; s2 < 3; ++s2) i.interactionStrength[s1][s2] = interaction.interactionStrength[s2][s1];
					}
					spinModel->interactions.push_back(std::pair<LatticeIterator, SpinModel::SpinInteraction>(lattice->fromParametrization(rid), i));
				}
			}
		}

		//print ldf file
		if (ldfPath != "")
		{
			//open file file
			std::ofstream ldfFile(ldfPath, std::ios::out);
			if (!ldfFile.is_open()) throw Exception(Exception::Type::IOError, "Could not write lattice debug information to file. ");
			ldfFile << "<lattice>" << std::endl;

			//write sites
			for (auto i1 = lattice->getRange(0); i1 != lattice->end(); ++i1)
			{
				geometry::Vec3<double> p = lattice->getSitePosition(i1);
				std::string parametrized = (i1 - lattice->begin() < lattice->size) ? "true" : "false";
				ldfFile << boost::format("\t<site id=\"%d\" x=\"%f\" y=\"%f\" z=\"%f\" parametrized=\"%s\"/>") % (i1 - lattice->begin()) % p.x % p.y % p.z % parametrized << std::endl;
			}

			//write bonds
			for (auto i1 = lattice->getRange(0); i1 != lattice->end(); ++i1)
			{
				for (auto i2 = lattice->getRange(0); i2 != lattice->end(); ++i2)
				{
					auto s1Parm = lattice->getSiteParameters(i1);
					auto s2Parm = lattice->getSiteParameters(i2);
					LatticeSite s1 = LatticeSite(std::get<0>(s1Parm), std::get<1>(s1Parm), std::get<2>(s1Parm), std::get<3>(s1Parm));
					LatticeSite s2 = LatticeSite(std::get<0>(s2Parm), std::get<1>(s2Parm), std::get<2>(s2Parm), std::get<3>(s2Parm));

					for (auto bond : uc.latticeBonds)
					{
						if (bond.isConnectingFromTo(s1, s2)) ldfFile << boost::format("\t<bond from=\"%d\" to=\"%d\" />") % (i1 - lattice->begin()) % (i2 - lattice->begin()) << std::endl;
					}
				}
			}

			//write interactions
			for (auto i1 = lattice->getRange(0); i1 != lattice->end(); ++i1)
			{
				for (auto i2 = lattice->getRange(0); i2 != lattice->end(); ++i2)
				{
					auto s1Parm = lattice->getSiteParameters(i1);
					auto s2Parm = lattice->getSiteParameters(i2);
					LatticeSite s1 = LatticeSite(std::get<0>(s1Parm), std::get<1>(s1Parm), std::get<2>(s1Parm), std::get<3>(s1Parm));
					LatticeSite s2 = LatticeSite(std::get<0>(s2Parm), std::get<1>(s2Parm), std::get<2>(s2Parm), std::get<3>(s2Parm));

					for (auto i : spinModelDefinition.interactions)
					{
						if (i.isConnectingFromTo(s1, s2))
						{
							ldfFile << boost::format("\t<interaction from=\"%d\" to=\"%d\" value=\"[[%f,%f,%f],[%f,%f,%f],[%f,%f,%f]]\" />") % (i1 - lattice->begin()) % (i2 - lattice->begin()) % i.interactionStrength[0][0] % i.interactionStrength[0][1] % i.interactionStrength[0][2] % i.interactionStrength[1][0] % i.interactionStrength[1][1] % i.interactionStrength[1][2] % i.interactionStrength[2][0] % i.interactionStrength[2][1] % i.interactionStrength[2][2] << std::endl;
						}
					}
				}
			}

			//finalize file
			ldfFile << "</lattice>" << std::endl;
			ldfFile.close();
		}

		//return lattice
		return std::pair<Lattice*, SpinModel*>(lattice, spinModel);
	}
}