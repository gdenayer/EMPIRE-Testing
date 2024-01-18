/******************************************************************************//**
 * \file main.cpp
 * This file holds the main function for testing mortar mapper.
 * \author Tianyang Wang
 * \date 9/21/2011
 *********************************************************************************/

#include <iostream>
#include <string>
#include "GiDFileIO.h"
#include "TriangleMortarMapper.h"
#include "FieldCreator.h"
#include "FluxCreator.h"
#include "QuadMortarMapper.h"
#include "AbstractMortarMapper.h"
#include "tinyxml.h"
#include <assert.h>

using namespace std;

double *copyDoubleArray(double *origin, int size);
int *copyIntArray(int *origin, int size);

static const int DOF_DIM = 3;

/*
 * ==================================================================================
 * =============================== main function ====================================
 * ==================================================================================
 */
int main(int argc, char** argv) {
	// 0. read input file by using tinyxml
	assert(argc == 2);
	char *file = argv[1];
	TiXmlDocument *doc = new TiXmlDocument(file);
	assert(doc->LoadFile());
	doc->Print();
	string masterFile = doc->RootElement()->FirstChild("meshes")->FirstChild(
			"master")->ToElement()->Attribute("path");
	string slaveFile = doc->RootElement()->FirstChild("meshes")->FirstChild(
			"slave")->ToElement()->Attribute("path");
	string mapperType =
			doc->RootElement()->FirstChild("mapper")->ToElement()->Attribute(
					"type");

	bool surfaceNormalOpposite = false;
	if (doc->RootElement()->FirstChild("geometry"))
		doc->RootElement()->FirstChild("geometry")->ToElement()->QueryBoolAttribute(
				"surfaceNormalOpposite", &surfaceNormalOpposite);
	bool triangleSymmetric = false;
	if (doc->RootElement()->FirstChild("mapper")->FirstChild("triangleMapper"))
		doc->RootElement()->FirstChild("mapper")->FirstChild("triangleMapper")->ToElement()->QueryBoolAttribute(
				"symmetric", &triangleSymmetric);

	string fieldFunc = "constant";
	string fluxFunc = "constant";
	bool perpendicularToSurface = false;
	if (doc->RootElement()->FirstChild("data")) {
		fieldFunc =
				doc->RootElement()->FirstChild("data")->FirstChild("field")->ToElement()->Attribute(
						"function");
		fluxFunc =
				doc->RootElement()->FirstChild("data")->FirstChild("flux")->ToElement()->Attribute(
						"function");
		doc->RootElement()->FirstChild("data")->FirstChild("flux")->ToElement()->QueryBoolAttribute(
				"perpendicularToSurface", &perpendicularToSurface);
	}

	delete doc;

	cout << "======================= Settings =====================" << endl;
	cout << "master mesh path         : " << masterFile << endl;
	cout << "slave mesh path          : " << slaveFile << endl;
	cout << "surface normal opposite  : " << surfaceNormalOpposite << endl;
	cout << "mapper                   : " << mapperType << endl;
	cout << "triangle mapper symmetric: " << triangleSymmetric << endl;
	cout << "field function           : " << fieldFunc << endl;
	cout << "flux function            : " << fluxFunc << endl;
	cout << "perpendicular to surface : " << perpendicularToSurface << endl;
	cout << "======================================================" << endl;

	// 1. extract data from MISlave
	int slaveNumNodes;
	int slaveNumElems;
	int slaveNodesPerElem;
	double *slaveNodeCoors;
	int *slaveNodeNumbers;
	int *slaveElemTable;
	GiDFileIO::readDotMsh(slaveFile, slaveNumNodes, slaveNumElems,
			slaveNodesPerElem, slaveNodeCoors, slaveNodeNumbers,
			slaveElemTable);

	// 2. extract data from MIMaster
	int masterNumNodes;
	int masterNumElems;
	int masterNodesPerElem;
	double *masterNodeCoors;
	int *masterNodeNumbers;
	int *masterElemTable;
	GiDFileIO::readDotMsh(masterFile, masterNumNodes, masterNumElems,
			masterNodesPerElem, masterNodeCoors, masterNodeNumbers,
			masterElemTable);

	// 3. generate DOFs to be mapped
	double *slaveFieldDOFs = new double[slaveNumNodes * DOF_DIM];
	double *masterFieldDOFs = new double[masterNumNodes * DOF_DIM];
	double *slaveFluxDOFs = new double[slaveNumNodes * DOF_DIM];
	double *masterFluxDOFs = new double[masterNumNodes * DOF_DIM];

	AbstractDataCreator *fieldSolver = new FieldCreator(slaveNumNodes, slaveNumElems,
			slaveNodesPerElem, slaveNodeCoors, slaveNodeNumbers, slaveElemTable,
			fieldFunc);
	fieldSolver->create(slaveFieldDOFs);

	AbstractDataCreator *fluxSolver = new FluxCreator(masterNumNodes, masterNumElems,
			masterNodesPerElem, masterNodeCoors, masterNodeNumbers,
			masterElemTable, fluxFunc, perpendicularToSurface);
	fluxSolver->create(masterFluxDOFs);

	// 4. do mortar consistent mapping
	// is it necessary to copy the coordinates before mapping?
	double *slaveNodeCoors_ = copyDoubleArray(slaveNodeCoors,
			slaveNumNodes * 3);
	int *slaveNodeNumbers_ = copyIntArray(slaveNodeNumbers, slaveNumNodes);
	int *slaveElemTable_ = copyIntArray(slaveElemTable,
			slaveNumElems * slaveNodesPerElem);
	double *masterNodeCoors_ = copyDoubleArray(masterNodeCoors,
			masterNumNodes * 3);
	int *masterNodeNumbers_ = copyIntArray(masterNodeNumbers, masterNumNodes);
	int *masterElemTable_ = copyIntArray(masterElemTable,
			masterNumElems * masterNodesPerElem);

	AbstractMortarMapper *mapper = NULL;
	if (mapperType == "quad")
		mapper = new QuadMortarMapper(slaveNumNodes, slaveNumElems,
				slaveNodesPerElem, slaveNodeCoors_, slaveNodeNumbers_,
				slaveElemTable_, masterNumNodes, masterNumElems,
				masterNodesPerElem, masterNodeCoors_, masterNodeNumbers_,
				masterElemTable_, surfaceNormalOpposite, false);
	else if (mapperType == "quad dual")
		mapper = new QuadMortarMapper(slaveNumNodes, slaveNumElems,
				slaveNodesPerElem, slaveNodeCoors_, slaveNodeNumbers_,
				slaveElemTable_, masterNumNodes, masterNumElems,
				masterNodesPerElem, masterNodeCoors_, masterNodeNumbers_,
				masterElemTable_, surfaceNormalOpposite, true);
	else if (mapperType == "triangle")
		mapper = new TriangleMortarMapper(slaveNumNodes, slaveNumElems,
				slaveNodesPerElem, slaveNodeCoors_, slaveNodeNumbers_,
				slaveElemTable_, masterNumNodes, masterNumElems,
				masterNodesPerElem, masterNodeCoors_, masterNodeNumbers_,
				masterElemTable_, surfaceNormalOpposite, false,
				triangleSymmetric);
	else if (mapperType == "triangle dual")
		mapper = new TriangleMortarMapper(slaveNumNodes, slaveNumElems,
				slaveNodesPerElem, slaveNodeCoors_, slaveNodeNumbers_,
				slaveElemTable_, masterNumNodes, masterNumElems,
				masterNodesPerElem, masterNodeCoors_, masterNodeNumbers_,
				masterElemTable_, surfaceNormalOpposite, true,
				triangleSymmetric);
	if (mapper == NULL) {
		cout << endl << "wrong mapper type!" << endl;
		exit(EXIT_FAILURE);
	}

	double *slaveFieldDOFs_x = new double[slaveNumNodes];
	double *slaveFieldDOFs_y = new double[slaveNumNodes];
	double *slaveFieldDOFs_z = new double[slaveNumNodes];
	for (int i = 0; i < slaveNumNodes; i++) {
		slaveFieldDOFs_x[i] = slaveFieldDOFs[i * DOF_DIM];
		slaveFieldDOFs_y[i] = slaveFieldDOFs[i * DOF_DIM + 1];
		slaveFieldDOFs_z[i] = slaveFieldDOFs[i * DOF_DIM + 2];
	}

	double *masterFieldDOFs_x = new double[masterNumNodes];
	double *masterFieldDOFs_y = new double[masterNumNodes];
	double *masterFieldDOFs_z = new double[masterNumNodes];

	mapper->consistentMapping(slaveFieldDOFs_x, masterFieldDOFs_x);
	mapper->consistentMapping(slaveFieldDOFs_y, masterFieldDOFs_y);
	mapper->consistentMapping(slaveFieldDOFs_z, masterFieldDOFs_z);

	for (int i = 0; i < masterNumNodes; i++) {
		masterFieldDOFs[i * DOF_DIM] = masterFieldDOFs_x[i];
		masterFieldDOFs[i * DOF_DIM + 1] = masterFieldDOFs_y[i];
		masterFieldDOFs[i * DOF_DIM + 2] = masterFieldDOFs_z[i];
	}

	// 5. do mortar conservative mapping
	double *masterFluxDOFs_x = new double[masterNumNodes];
	double *masterFluxDOFs_y = new double[masterNumNodes];
	double *masterFluxDOFs_z = new double[masterNumNodes];
	for (int i = 0; i < masterNumNodes; i++) {
		masterFluxDOFs_x[i] = masterFluxDOFs[i * DOF_DIM];
		masterFluxDOFs_y[i] = masterFluxDOFs[i * DOF_DIM + 1];
		masterFluxDOFs_z[i] = masterFluxDOFs[i * DOF_DIM + 2];
	}

	double *slaveFluxDOFs_x = new double[slaveNumNodes];
	double *slaveFluxDOFs_y = new double[slaveNumNodes];
	double *slaveFluxDOFs_z = new double[slaveNumNodes];
	mapper->conservativeMapping(masterFluxDOFs_x, slaveFluxDOFs_x);
	mapper->conservativeMapping(masterFluxDOFs_y, slaveFluxDOFs_y);
	mapper->conservativeMapping(masterFluxDOFs_z, slaveFluxDOFs_z);

	for (int i = 0; i < slaveNumNodes; i++) {
		slaveFluxDOFs[i * DOF_DIM] = slaveFluxDOFs_x[i];
		slaveFluxDOFs[i * DOF_DIM + 1] = slaveFluxDOFs_y[i];
		slaveFluxDOFs[i * DOF_DIM + 2] = slaveFluxDOFs_z[i];
	}

	// 6. Error computation
	double *masterError = new double[masterNumNodes * DOF_DIM];

	AbstractDataCreator *fieldError = new FieldCreator(masterNumNodes, masterNumElems,
			masterNodesPerElem, masterNodeCoors, masterNodeNumbers,
			masterElemTable, fieldFunc);
	fieldError->create(masterError);
	for (int i = 0; i < masterNumNodes * DOF_DIM; i++)
		masterError[i] = fabs(masterError[i] - masterFieldDOFs[i]);

	double *slaveError = new double[slaveNumNodes * DOF_DIM];

	AbstractDataCreator *fluxError = new FluxCreator(slaveNumNodes, slaveNumElems,
			slaveNodesPerElem, slaveNodeCoors, slaveNodeNumbers, slaveElemTable,
			fluxFunc, perpendicularToSurface);
	fluxError->create(slaveError);
	for (int i = 0; i < slaveNumNodes * DOF_DIM; i++)
		slaveError[i] = fabs(slaveError[i] - slaveFluxDOFs[i]);

	// 7. write output
	string analysisName = "\"Mortar_mapping\"";
	const int NUM_RESULTS = 3;

	string slaveMesh = slaveFile;
	slaveMesh.erase(slaveMesh.end() - 4, slaveMesh.end());
	string slaveResultFile = slaveMesh;
	slaveResultFile.append(".res");
	string slaveResultNames[NUM_RESULTS] = { "\"slave field\"",
			"\"slave flux\"", "\"slave flux error\"" };
	double *slaveResults[NUM_RESULTS] = { slaveFieldDOFs, slaveFluxDOFs,
			slaveError };
	GiDFileIO::writeDotRes(slaveResultFile, slaveNumNodes, slaveNodeNumbers,
			analysisName, NUM_RESULTS, slaveResultNames, slaveResults);
	GiDFileIO::writeGnuplotError(slaveMesh, slaveNumNodes, slaveNodeCoors,
			slaveError);

	string masterMesh = masterFile;
	masterMesh.erase(masterMesh.end() - 4, masterMesh.end());
	string masterResultFile = masterMesh;
	masterResultFile.append(".res");
	string masterResultNames[NUM_RESULTS] = { "\"master field\"",
			"\"master flux\"", "\"master field error\"" };
	double *masterResults[NUM_RESULTS] = { masterFieldDOFs, masterFluxDOFs,
			masterError };
	GiDFileIO::writeDotRes(masterResultFile, masterNumNodes, masterNodeNumbers,
			analysisName, NUM_RESULTS, masterResultNames, masterResults);
	GiDFileIO::writeGnuplotError(masterMesh, masterNumNodes, masterNodeCoors,
			masterError);

	//delete all
	delete mapper;
	cout << endl << endl << "Bye!" << endl;
	return 0;
}

/*
 * ==================================================================================
 * ========================== auxiliary functions ===================================
 * ==================================================================================
 */

double *copyDoubleArray(double *origin, int size) {
	double *copy = new double[size];
	for (int i = 0; i < size; i++)
		copy[i] = origin[i];
	return copy;
}

int *copyIntArray(int *origin, int size) {
	int *copy = new int[size];
	for (int i = 0; i < size; i++)
		copy[i] = origin[i];
	return copy;
}
