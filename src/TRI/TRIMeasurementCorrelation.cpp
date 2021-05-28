/**
 * @file TRIMeasurementCorrelation.cpp
 * @author Finn Lasse Buessen
 * @brief Correlation measurement for time reversal invariant models.
 * 
 * @copyright Copyright (c) 2020
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include <hdf5.h>
#include "lib/Integrator.hpp"
#include "lib/ValueBundle.hpp"
#include "TRIMeasurementCorrelation.hpp"
#include "SpinParser.hpp"
#include "TRIFrgCore.hpp"
#include "TRIEffectiveAction.hpp"

TRIMeasurementCorrelation::TRIMeasurementCorrelation(const std::string &outfile, const float minCutoff, const float maxCutoff, const bool defer) : Measurement(outfile, minCutoff, maxCutoff, defer, true)
{
	_currentCutoff = -1.0f;
	int latticeSizeExtended = 0;
	for (auto i = FrgCommon::lattice().getRange(0); i != FrgCommon::lattice().end(); ++i) ++latticeSizeExtended;
	int latticeSizeBasis = int(FrgCommon::lattice()._basis.size());
	_memoryStepLattice = latticeSizeBasis * latticeSizeExtended;

	//prepare correlation buffer
	_correlationsDD = new float[latticeSizeBasis * latticeSizeExtended];
	_correlationsXX = new float[latticeSizeBasis * latticeSizeExtended];
	_correlationsXY = new float[latticeSizeBasis * latticeSizeExtended];
	_correlationsXZ = new float[latticeSizeBasis * latticeSizeExtended];
	_correlationsYX = new float[latticeSizeBasis * latticeSizeExtended];
	_correlationsYY = new float[latticeSizeBasis * latticeSizeExtended];
	_correlationsYZ = new float[latticeSizeBasis * latticeSizeExtended];
	_correlationsZX = new float[latticeSizeBasis * latticeSizeExtended];
	_correlationsZY = new float[latticeSizeBasis * latticeSizeExtended];
	_correlationsZZ = new float[latticeSizeBasis * latticeSizeExtended];

	//set up loadManager
	//stack0
	HMP::StackIdentifier dataStack0 = SpinParser::spinParser()->getLoadManager()->addMasterStackExplicit<float>(
		&_currentCutoff,
		1,
		[](int linearIterator)->float { return SpinParser::spinParser()->getFrgCore()->flowingFunctional()->cutoff; },
		1,
		1,
		true);
	//stack1
	HMP::StackIdentifier dataStack1 = SpinParser::spinParser()->getLoadManager()->addMasterStackImplicit<float>(
		_correlationsXX,
		1,
		std::bind(&TRIMeasurementCorrelation::_calculateCorrelation, this, std::placeholders::_1),
		latticeSizeBasis * latticeSizeExtended,
		1,
		1,
		false);
	//stack2
	SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		_correlationsXY,
		1,
		dataStack1,
		latticeSizeBasis * latticeSizeExtended);
	//stack3
	SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		_correlationsXZ,
		1,
		dataStack1,
		latticeSizeBasis * latticeSizeExtended);
	//stack4
	SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		_correlationsYX,
		1,
		dataStack1,
		latticeSizeBasis * latticeSizeExtended);
	//stack5
	SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		_correlationsYY,
		1,
		dataStack1,
		latticeSizeBasis * latticeSizeExtended);
	//stack6
	SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		_correlationsYZ,
		1,
		dataStack1,
		latticeSizeBasis * latticeSizeExtended);
	//stack7
	SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		_correlationsZX,
		1,
		dataStack1,
		latticeSizeBasis * latticeSizeExtended);
	//stack8
	SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		_correlationsZY,
		1,
		dataStack1,
		latticeSizeBasis * latticeSizeExtended);
	//stack9
	SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		_correlationsZZ,
		1,
		dataStack1,
		latticeSizeBasis * latticeSizeExtended);
	//stack10
	SpinParser::spinParser()->getLoadManager()->addSlaveStack<float>(
		_correlationsDD,
		1,
		dataStack1,
		latticeSizeBasis * latticeSizeExtended);
	_loadManagedStacks.insert(_loadManagedStacks.end(), { dataStack0, dataStack1 });
}

TRIMeasurementCorrelation::~TRIMeasurementCorrelation()
{
	delete[] _correlationsXX;
	delete[] _correlationsXY;
	delete[] _correlationsXZ;
	delete[] _correlationsYX;
	delete[] _correlationsYY;
	delete[] _correlationsYZ;
	delete[] _correlationsZX;
	delete[] _correlationsZY;
	delete[] _correlationsZZ;
	delete[] _correlationsDD;
}

void TRIMeasurementCorrelation::takeMeasurement(const EffectiveAction &state, const bool isMasterTask) const
{
	if (_currentCutoff != state.cutoff) SpinParser::spinParser()->getLoadManager()->calculate(_loadManagedStacks.data(), int(_loadManagedStacks.size()));

	if (isMasterTask)
	{
		_writeOutfileCorrelation("TRICorXX", _correlationsXX);
		_writeOutfileCorrelation("TRICorXY", _correlationsXY);
		_writeOutfileCorrelation("TRICorXZ", _correlationsXZ);
		_writeOutfileCorrelation("TRICorYX", _correlationsYX);
		_writeOutfileCorrelation("TRICorYY", _correlationsYY);
		_writeOutfileCorrelation("TRICorYZ", _correlationsYZ);
		_writeOutfileCorrelation("TRICorZX", _correlationsZX);
		_writeOutfileCorrelation("TRICorZY", _correlationsZY);
		_writeOutfileCorrelation("TRICorZZ", _correlationsZZ);
		_writeOutfileCorrelation("TRICorDD", _correlationsDD);
	}
}

void TRIMeasurementCorrelation::_calculateCorrelation(const int iterator) const
{
	//calculate real space susceptibility
	float nu = 0.0f;
	float cut = SpinParser::spinParser()->getFrgCore()->flowingFunctional()->cutoff;
	TRIVertexSingleParticle *v2 = static_cast<TRIEffectiveAction *>(SpinParser::spinParser()->getFrgCore()->flowingFunctional())->vertexSingleParticle;
	TRIVertexTwoParticle *v4 = static_cast<TRIEffectiveAction *>(SpinParser::spinParser()->getFrgCore()->flowingFunctional())->vertexTwoParticle;

	ValueSuperbundle<float, 16> susceptibility(FrgCommon::lattice().size);
	ValueSuperbundle<float, 16> stackBuffer(FrgCommon::lattice().size);
	ValueSuperbundle<float, 16> buffer1(FrgCommon::lattice().size);
	ValueSuperbundle<float, 16> buffer2(FrgCommon::lattice().size);
	ValueSuperbundle<float, 16> buffer3(FrgCommon::lattice().size);
	ValueSuperbundle<float, 16> buffer4(FrgCommon::lattice().size);

	//integration kernel
	std::function<void(float, ValueSuperbundle<float, 16> &)> integralKernel = [&](const float w, ValueSuperbundle<float, 16> &returnBuffer) -> void
	{
		returnBuffer.reset();

		//term1
		float term1 = 1.0f / ((w + v2->getValue(w)) * (w + nu + v2->getValue(w + nu)));
		returnBuffer.bundle(15)[0] += 2.0f * term1 / float(2.0f * M_PI);
		returnBuffer.bundle(0)[0] += 0.5f * term1 / float(2.0f * M_PI);
		returnBuffer.bundle(5)[0] += 0.5f * term1 / float(2.0f * M_PI);
		returnBuffer.bundle(10)[0] += 0.5f * term1 / float(2.0f * M_PI);

		//term2
		std::function<void(float, ValueSuperbundle<float, 16> &)> innerKernel = [&](const float wp, ValueSuperbundle<float, 16> &ret) -> void
		{
			ret.reset();
			const TRIVertexTwoParticleAccessBuffer<8> ab0 = v4->generateAccessBuffer(w + wp + nu, nu, w - wp);
			v4->getValueSuperbundle(ab0, stackBuffer);

			const TRIVertexTwoParticleAccessBuffer<8> ab1 = v4->generateAccessBuffer(w + wp + nu, w - wp, nu);
			const float vxx = v4->getValueLocal(SpinComponent::X, SpinComponent::X, ab1);
			const float vxy = v4->getValueLocal(SpinComponent::X, SpinComponent::Y, ab1);
			const float vxz = v4->getValueLocal(SpinComponent::X, SpinComponent::Z, ab1);
			const float vxd = v4->getValueLocal(SpinComponent::X, SpinComponent::None, ab1);
			const float vyx = v4->getValueLocal(SpinComponent::Y, SpinComponent::X, ab1);
			const float vyy = v4->getValueLocal(SpinComponent::Y, SpinComponent::Y, ab1);
			const float vyz = v4->getValueLocal(SpinComponent::Y, SpinComponent::Z, ab1);
			const float vyd = v4->getValueLocal(SpinComponent::Y, SpinComponent::None, ab1);
			const float vzx = v4->getValueLocal(SpinComponent::Z, SpinComponent::X, ab1);
			const float vzy = v4->getValueLocal(SpinComponent::Z, SpinComponent::Y, ab1);
			const float vzz = v4->getValueLocal(SpinComponent::Z, SpinComponent::Z, ab1);
			const float vzd = v4->getValueLocal(SpinComponent::Z, SpinComponent::None, ab1);
			const float vdx = v4->getValueLocal(SpinComponent::None, SpinComponent::X, ab1);
			const float vdy = v4->getValueLocal(SpinComponent::None, SpinComponent::Y, ab1);
			const float vdz = v4->getValueLocal(SpinComponent::None, SpinComponent::Z, ab1);
			const float vdd = v4->getValueLocal(SpinComponent::None, SpinComponent::None, ab1);

			//dumbbell diagram
			ret.bundle(15).multSub(4.0f, stackBuffer.bundle(15));
			ret.bundle(0).multSub(1.0f, stackBuffer.bundle(0));
			ret.bundle(1).multSub(1.0f, stackBuffer.bundle(1));
			ret.bundle(2).multSub(1.0f, stackBuffer.bundle(2));
			ret.bundle(4).multSub(1.0f, stackBuffer.bundle(4));
			ret.bundle(5).multSub(1.0f, stackBuffer.bundle(5));
			ret.bundle(6).multSub(1.0f, stackBuffer.bundle(6));
			ret.bundle(8).multSub(1.0f, stackBuffer.bundle(8));
			ret.bundle(9).multSub(1.0f, stackBuffer.bundle(9));
			ret.bundle(10).multSub(1.0f, stackBuffer.bundle(10));

			//egg diagram
			ret.bundle(15)[0] += 2.0f * vdd;
			ret.bundle(15)[0] += 2.0f * vzz;
			ret.bundle(15)[0] += 2.0f * vyy;
			ret.bundle(15)[0] += 2.0f * vxx;
			ret.bundle(0)[0] += 0.5f * vdd;
			ret.bundle(0)[0] -= 0.5f * vzz;
			ret.bundle(0)[0] -= 0.5f * vyy;
			ret.bundle(0)[0] += 0.5f * vxx;
			ret.bundle(1)[0] += 0.5f * vdz;
			ret.bundle(1)[0] -= 0.5f * vzd;
			ret.bundle(1)[0] += 0.5f * vyx;
			ret.bundle(1)[0] += 0.5f * vxy;
			ret.bundle(2)[0] -= 0.5f * vdy;
			ret.bundle(2)[0] += 0.5f * vzx;
			ret.bundle(2)[0] += 0.5f * vyd;
			ret.bundle(2)[0] += 0.5f * vxz;
			ret.bundle(4)[0] -= 0.5f * vdz;
			ret.bundle(4)[0] += 0.5f * vzd;
			ret.bundle(4)[0] += 0.5f * vyx;
			ret.bundle(4)[0] += 0.5f * vxy;
			ret.bundle(5)[0] += 0.5f * vdd;
			ret.bundle(5)[0] -= 0.5f * vzz;
			ret.bundle(5)[0] += 0.5f * vyy;
			ret.bundle(5)[0] -= 0.5f * vxx;
			ret.bundle(6)[0] += 0.5f * vdx;
			ret.bundle(6)[0] += 0.5f * vzy;
			ret.bundle(6)[0] += 0.5f * vyz;
			ret.bundle(6)[0] -= 0.5f * vxd;
			ret.bundle(8)[0] += 0.5f * vdy;
			ret.bundle(8)[0] += 0.5f * vzx;
			ret.bundle(8)[0] -= 0.5f * vyd;
			ret.bundle(8)[0] += 0.5f * vxz;
			ret.bundle(9)[0] -= 0.5f * vdx;
			ret.bundle(9)[0] += 0.5f * vzy;
			ret.bundle(9)[0] += 0.5f * vyz;
			ret.bundle(9)[0] += 0.5f * vxd;
			ret.bundle(10)[0] += 0.5f * vdd;
			ret.bundle(10)[0] += 0.5f * vzz;
			ret.bundle(10)[0] -= 0.5f * vyy;
			ret.bundle(10)[0] -= 0.5f * vxx;

			float normalization = 1.0f / ((w + v2->getValue(w)) * (w + nu + v2->getValue(w + nu)) * (wp + v2->getValue(wp)) * (wp + nu + v2->getValue(wp + nu)) * float(4.0f * M_PI * M_PI));
			ret *= normalization;
		};
		if (-(nu + cut) > *FrgCommon::frequency().beginNegative())
		{
			ImplicitIntegrator::integrateWithObscureRightBoundary(FrgCommon::frequency().beginNegative(), -cut - nu, innerKernel, buffer3, buffer4);
			returnBuffer += buffer4;
		}
		if (nu - cut > cut)
		{
			ImplicitIntegrator::integrateWithObscureBoundaries(-nu + cut, -cut, innerKernel, buffer3, buffer4);
			returnBuffer += buffer4;
		}
		if (cut < *FrgCommon::frequency().last())
		{
			ImplicitIntegrator::integrateWithObscureLeftBoundary(cut, FrgCommon::frequency().last(), innerKernel, buffer3, buffer4);
			returnBuffer += buffer4;
		}
	};

	if (-(nu + cut) > *FrgCommon::frequency().beginNegative())
	{
		ImplicitIntegrator::integrateWithObscureRightBoundary(FrgCommon::frequency().beginNegative(), -cut - nu, integralKernel, buffer1, buffer2);
		susceptibility += buffer2;
	}
	if (nu - cut > cut)
	{
		ImplicitIntegrator::integrateWithObscureBoundaries(-nu + cut, -cut, integralKernel, buffer1, buffer2);
		susceptibility += buffer2;
	}
	if (cut < *FrgCommon::frequency().last())
	{
		ImplicitIntegrator::integrateWithObscureLeftBoundary(cut, FrgCommon::frequency().last(), integralKernel, buffer1, buffer2);
		susceptibility += buffer2;
	}

	int offset = iterator * _memoryStepLattice;
	for (auto i = FrgCommon::lattice().getBasis(); i != FrgCommon::lattice().end(); ++i)
	{
		for (auto j = FrgCommon::lattice().getRange(i); j != FrgCommon::lattice().end(); ++j)
		{
			int rid;

			SpinComponent sx(SpinComponent::X);
			SpinComponent sy(SpinComponent::Y);
			SpinComponent sz(SpinComponent::Z);
			rid = FrgCommon::lattice().symmetryTransform(i, j, sx, sy, sz);

			_correlationsDD[offset] = susceptibility.bundle(15)[rid];
			_correlationsXX[offset] = susceptibility.bundle(4 * static_cast<int>(sx) + static_cast<int>(sx))[rid];
			_correlationsXY[offset] = susceptibility.bundle(4 * static_cast<int>(sx) + static_cast<int>(sy))[rid];
			_correlationsXZ[offset] = susceptibility.bundle(4 * static_cast<int>(sx) + static_cast<int>(sz))[rid];
			_correlationsYX[offset] = susceptibility.bundle(4 * static_cast<int>(sy) + static_cast<int>(sx))[rid];
			_correlationsYY[offset] = susceptibility.bundle(4 * static_cast<int>(sy) + static_cast<int>(sy))[rid];
			_correlationsYZ[offset] = susceptibility.bundle(4 * static_cast<int>(sy) + static_cast<int>(sz))[rid];
			_correlationsZX[offset] = susceptibility.bundle(4 * static_cast<int>(sz) + static_cast<int>(sx))[rid];
			_correlationsZY[offset] = susceptibility.bundle(4 * static_cast<int>(sz) + static_cast<int>(sy))[rid];
			_correlationsZZ[offset] = susceptibility.bundle(4 * static_cast<int>(sz) + static_cast<int>(sz))[rid];

			++offset;
		}
	}
}

void TRIMeasurementCorrelation::_writeOutfileHeader(const std::string &observableGroup) const
{
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	//open file
	hid_t file = H5Fopen(outfile().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	if (file < 0) throw Exception(Exception::Type::IOError, "Could not open observable file [" + outfile() + "] for writing.");

	//open group
	hid_t group = H5Gopen(file, observableGroup.c_str(), H5P_DEFAULT);
	if (group < 0) throw Exception(Exception::Type::IOError, "Could not open obsfile group [" + observableGroup + "] for writing. ");

	//create meta group
	if (H5Lexists(group, "meta", H5P_DEFAULT) == 0)
	{
		hid_t mgroup = H5Gcreate(group, "meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		const int dataTypeLatticeSiteDim = 1;
		const hsize_t dataTypeLatticeSiteSize[dataTypeLatticeSiteDim] = { 3 };
		hid_t dataTypeLatticeSite = H5Tarray_create(H5T_NATIVE_FLOAT, dataTypeLatticeSiteDim, dataTypeLatticeSiteSize);

		//write lattice vectors
		float *latticeBuffer = new float[3 * 3];
		int i = 0;
		for (auto a = FrgCommon::lattice()._bravaisLattice.begin(); a != FrgCommon::lattice()._bravaisLattice.end(); ++a)
		{
			latticeBuffer[3 * i] = float(a->x);
			latticeBuffer[3 * i + 1] = float(a->y);
			latticeBuffer[3 * i + 2] = float(a->z);
			++i;
		}
		const int attrSpaceDimLattice = 1;
		const hsize_t attrSpaceSizeLattice[attrSpaceDimLattice] = { FrgCommon::lattice()._bravaisLattice.size() };
		hid_t attrSpaceLattice = H5Screate_simple(attrSpaceDimLattice, attrSpaceSizeLattice, NULL);
		hid_t datasetLattice = H5Dcreate(mgroup, "latticeVectors", dataTypeLatticeSite, attrSpaceLattice, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(datasetLattice, dataTypeLatticeSite, H5S_ALL, H5S_ALL, H5P_DEFAULT, latticeBuffer);
		H5Dclose(datasetLattice);
		H5Sclose(attrSpaceLattice);
		delete[] latticeBuffer;

		//write basis
		float *basisBuffer = new float[3 * FrgCommon::lattice()._basis.size()];
		i = 0;
		for (auto b = FrgCommon::lattice().getBasis(); b != FrgCommon::lattice().end(); ++b)
		{
			basisBuffer[3 * i] = float(FrgCommon::lattice().getSitePosition(b).x);
			basisBuffer[3 * i + 1] = float(FrgCommon::lattice().getSitePosition(b).y);
			basisBuffer[3 * i + 2] = float(FrgCommon::lattice().getSitePosition(b).z);
			++i;
		}
		const int attrSpaceDimBasis = 1;
		const hsize_t attrSpaceSizeBasis[attrSpaceDimBasis] = { FrgCommon::lattice()._basis.size() };
		hid_t attrSpaceBasis = H5Screate_simple(attrSpaceDimBasis, attrSpaceSizeBasis, NULL);
		hid_t datasetBasis = H5Dcreate(mgroup, "basis", dataTypeLatticeSite, attrSpaceBasis, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(datasetBasis, dataTypeLatticeSite, H5S_ALL, H5S_ALL, H5P_DEFAULT, basisBuffer);
		H5Dclose(datasetBasis);
		H5Sclose(attrSpaceBasis);
		delete[] basisBuffer;

		//write sites
		int inRangeCount = 0;
		for (SublatticeIterator i = FrgCommon::lattice().getRange(0); i != FrgCommon::lattice().end(); ++i) ++inRangeCount;

		const int dataSpaceDimSitesReference = 2;
		const hsize_t dataSpaceSizeSitesReference[dataSpaceDimSitesReference] = { FrgCommon::lattice()._basis.size(), hsize_t(inRangeCount) };
		hid_t dataSpaceSitesReference = H5Screate_simple(dataSpaceDimSitesReference, dataSpaceSizeSitesReference, NULL);

		float *SitesReferenceBuffer = new float[FrgCommon::lattice()._basis.size() * inRangeCount * 3];
		int j = 0;
		for (unsigned int b = 0; b < FrgCommon::lattice()._basis.size(); ++b)
		{
			for (SublatticeIterator i = FrgCommon::lattice().getRange(b); i != FrgCommon::lattice().end(); ++i)
			{
				*(SitesReferenceBuffer + 3 * j) = (float)FrgCommon::lattice().getSitePosition(i).x;
				*(SitesReferenceBuffer + 3 * j + 1) = (float)FrgCommon::lattice().getSitePosition(i).y;
				*(SitesReferenceBuffer + 3 * j + 2) = (float)FrgCommon::lattice().getSitePosition(i).z;
				++j;
			}
		}
		hid_t datasetSitesReference = H5Dcreate(mgroup, "sites", dataTypeLatticeSite, dataSpaceSitesReference, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(datasetSitesReference, dataTypeLatticeSite, H5S_ALL, H5S_ALL, H5P_DEFAULT, SitesReferenceBuffer);
		H5Dclose(datasetSitesReference);
		H5Sclose(dataSpaceSitesReference);
		delete[] SitesReferenceBuffer;

		//close meta group
		H5Tclose(dataTypeLatticeSite);
		H5Gclose(mgroup);
	}
	else
	{
		Log::log << Log::LogLevel::Warning << "The observable output file [" + outfile() + "] already contains the group [" + observableGroup + "/meta]. Skipping writing this information. " << Log::endl;
	}


	//close file
	H5Gclose(group);
	H5Fclose(file);
}

void TRIMeasurementCorrelation::_writeOutfileCorrelation(const std::string &observableGroup, const float *correlation) const
{
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	//open or create file
	hid_t file = (H5Fis_hdf5(outfile().c_str()) > 0) ? H5Fopen(outfile().c_str(), H5F_ACC_RDWR, H5P_DEFAULT) : H5Fcreate(outfile().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file < 0) throw Exception(Exception::Type::IOError, "Could not open observable file [" + outfile() + "] for writing");

	//open or create group
	hid_t group = (H5Lexists(file, observableGroup.c_str(), H5P_DEFAULT) > 0) ? H5Gopen(file, observableGroup.c_str(), H5P_DEFAULT) : H5Gcreate(file, observableGroup.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (group < 0) throw Exception(Exception::Type::IOError, "Could not open obsfile group [" + observableGroup + "] for writing");

	//ensure that meta information is included
	if (H5Lexists(group, "meta", H5P_DEFAULT) == 0) _writeOutfileHeader(observableGroup);

	//open or create data collection
	hid_t data = (H5Lexists(group, "data", H5P_DEFAULT) > 0) ? H5Gopen(group, "data", H5P_DEFAULT) : H5Gcreate(group, "data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (data < 0) throw Exception(Exception::Type::IOError, "Could not open obsfile group [" + observableGroup + "/data] for writing");

	//determine unique dataset name and check for duplicate dataset
	int datasetId = 0;
	hsize_t numDatasets;
	H5Gget_num_objs(data, &numDatasets);
	for (int i = 0; i < int(numDatasets); ++i)
	{
		if (H5Gget_objtype_by_idx(data, i) == H5G_GROUP)
		{
			++datasetId;

			const int datasetNameMaxLength = 32;
			char datasetName[datasetNameMaxLength];
			H5Gget_objname_by_idx(data, i, datasetName, datasetNameMaxLength);

			hid_t measurement = H5Gopen(data, datasetName, H5P_DEFAULT);
			hid_t attr = H5Aopen(measurement, "cutoff", H5P_DEFAULT);
			float c;
			H5Aread(attr, H5T_NATIVE_FLOAT, &c);
			H5Aclose(attr);
			H5Gclose(measurement);

			if (c == _currentCutoff)
			{
				Log::log << Log::LogLevel::Warning << "Found existing correlation measurement at cutoff " + std::to_string(_currentCutoff) + ". Discarding duplicate entry." << Log::endl;
				return;
			}
		}
	}
	std::string datasetName = "measurement_" + std::to_string(datasetId);

	//create new measurement group
	hid_t measurement = H5Gcreate(data, datasetName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	const int attrSpaceDim = 1;
	const hsize_t attrSpaceSize[1] = { 1 };
	hid_t attrSpace = H5Screate_simple(attrSpaceDim, attrSpaceSize, NULL);
	hid_t attr = H5Acreate(measurement, "cutoff", H5T_NATIVE_FLOAT, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_NATIVE_FLOAT, &_currentCutoff);
	H5Aclose(attr);
	H5Sclose(attrSpace);

	//write data
	int inRangeCount = 0;
	for (SublatticeIterator i = FrgCommon::lattice().getRange(0); i != FrgCommon::lattice().end(); ++i) ++inRangeCount;
	const int dataSpaceDim = 2;
	const hsize_t dataSpaceSize[2] = { FrgCommon::lattice()._basis.size(), hsize_t(inRangeCount) };
	hid_t dataSpace = H5Screate_simple(dataSpaceDim, dataSpaceSize, NULL);

	hid_t dataset = H5Dcreate(measurement, "data", H5T_NATIVE_FLOAT, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, correlation);
	H5Dclose(dataset);

	H5Sclose(dataSpace);

	//clean up
	H5Gclose(measurement);
	H5Gclose(data);
	H5Gclose(group);
	H5Fclose(file);
}