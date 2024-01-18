#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <assert.h>
#include "GiDMeshInventory.h"
#include "GiDBasicData.h"

GiDMeshInventory::GiDMeshInventory(vector<GiDMesh*> *vecMesh) :
		_vecMesh(vecMesh) {
}

GiDMeshInventory::GiDMeshInventory(string fileName) {
	ifstream meshFile(fileName.c_str());
	if (!meshFile) {
		cerr << "GiDMeshInventory::GiDMeshInventory: file \"" << fileName
				<< "\" cannot be found" << '\n';
		exit(EXIT_FAILURE);
	}

	/* parse the header line ***********************************************************************/
	string textLine;
	string dummy;
	string meshName;
	string elementType;
	int dimension;
	int nodesPerElem;
	istringstream *lineStream = new istringstream("");

	_vecMesh = new vector<GiDMesh*>();

	/* 1) ALL meshes must define their own coordinates, forbid shared coordinates
	 * 2) No comment line between 'coordinates' and 'end coordinates', or 'elements' and 'end elements'
	 */
	while (getline(meshFile, textLine, '\n')) {
		delete lineStream;
		lineStream = new istringstream(textLine.c_str());
		(*lineStream) >> dummy;

		if (strToLower(dummy).compare("mesh") == 0) {

			int count = 0;
			while ((*lineStream) >> dummy)
				count++;

			const int headerLineWordsWithoutMeshName = 6;
			if (count > headerLineWordsWithoutMeshName) { // there is a mesh name
				delete lineStream;
				lineStream = new istringstream(textLine.c_str());
				(*lineStream) >> dummy;
				do {
					(*lineStream) >> dummy;
					meshName.append(" ").append(dummy);
				} while (meshName[meshName.size() - 1] != '\"'); // parse the name with spaces in
				meshName.erase(meshName.begin());

				(*lineStream) >> dummy; // skip the keyword "dimension"

				(*lineStream) >> dimension >> dummy >> elementType >> dummy
						>> nodesPerElem;
			} else { // there is no mesh name
				delete lineStream;
				lineStream = new istringstream(textLine.c_str());
				(*lineStream) >> dummy;
				(*lineStream) >> dummy;
				(*lineStream) >> dimension >> dummy >> elementType >> dummy
						>> nodesPerElem;
			}

			/* parse the coordinates paragraph **************************************************************************/
			do {
				getline(meshFile, textLine, '\n');
				delete lineStream;
				lineStream = new istringstream(textLine.c_str());
				(*lineStream) >> dummy;
			} while (strToLower(dummy).compare("coordinates") != 0); // find the line starting with keyword "coordinates"

			char e;
			int id;
			vector<GiDNode*> *vecNode = new vector<GiDNode*>();
			meshFile >> e;
			while (tolower(e) != 'e') {
				meshFile.putback(e);
				meshFile >> id;
				double *coor = new double[dimension];
				for (int i = 0; i < dimension; i++) {
					meshFile >> coor[i];
				}
				GiDNode *node = new GiDNode(id, coor);
				vecNode->push_back(node);
				meshFile >> e;
			}

			/* parse the elements paragraph **********************************************************************/
			do {
				getline(meshFile, textLine, '\n');
				delete lineStream;
				lineStream = new istringstream(textLine.c_str());
				(*lineStream) >> dummy;
			} while (strToLower(dummy).compare("elements") != 0); // find the line starting with keyword "elements"

			bool hasMat;
			streampos curPos = meshFile.tellg();
			// this block judges whether the elements have a material number
			{
				getline(meshFile, textLine, '\n');
				delete lineStream;
				lineStream = new istringstream(textLine.c_str());
				int tmp;
				int i = 0;
				while ((*lineStream) >> tmp) {
					i++;
				}
				if (i == (nodesPerElem + 1))
					hasMat = false;
				else
					hasMat = true;
			}
			assert(!hasMat);
			meshFile.seekg(curPos);

			vector<GiDElement*> *vecElem = new vector<GiDElement*>();

			meshFile >> e;
			while (tolower(e) != 'e') {
				meshFile.putback(e);
				meshFile >> id;
				int *nodeIds = new int[nodesPerElem];
				for (int i = 0; i < nodesPerElem; i++) {
					meshFile >> nodeIds[i];
				}
				GiDElement *element = new GiDElement(id, nodeIds);
				vecElem->push_back(element);
				meshFile >> e;
			}

			GiDMesh* mesh = new GiDMesh(meshName, elementType, dimension,
					nodesPerElem, vecNode, vecElem);
			_vecMesh->push_back(mesh);
		}
	}
}

GiDMeshInventory::~GiDMeshInventory() {
	for (unsigned i = 0; i < _vecMesh->size(); i++)
		delete _vecMesh->at(i);
	delete _vecMesh;
}

vector<GiDMesh*> *GiDMeshInventory::getMeshes() {
	return _vecMesh;
}

void GiDMeshInventory::writeDotMsh(ostream &out) {
	out << "#Edited by GiDParser:\n\n";
	vector<GiDMesh*> *meshes = getMeshes();
	for (unsigned j = 0; j < (*meshes).size(); j++) {
		GiDMesh *mesh = (*meshes)[j];
		out << "mesh " << mesh->getName() << " dimension "
				<< mesh->getDimension() << " elemType "
				<< mesh->getElementType() << " Nnode "
				<< mesh->getNodesPerElem() << '\n';
		out << "coordinates" << '\n';
		vector<GiDNode*> *vecNode = mesh->getNodes();
		int numNodes = vecNode->size();
		GiDNode *node;
		for (int i = 0; i < numNodes; i++) {
			node = (*vecNode)[i];
			out << setw(5) << node->getId() << "   ";
			for (int j = 0; j < mesh->getDimension(); j++)
				out << setw(15) << (node->getCoordinate())[j];
			out << '\n';
		}
		out << "end coordinates\n\n";

		out << "elements" << '\n';
		vector<GiDElement*> *vecElem = mesh->getElements();
		int numElems = vecElem->size();
		GiDElement *elem;
		for (int i = 0; i < numElems; i++) {
			elem = (*vecElem)[i];
			out << setw(5) << elem->getId() << "   ";
			for (int j = 0; j < mesh->getNodesPerElem(); j++)
				out << setw(8) << (elem->getNodeIds())[j];
			out << '\n';
		}
		out << "end elements\n\n";
	}
}

string GiDMeshInventory::strToLower(string strInput) {
	string strLower;
	int i = 0;
	while (strInput[i]) {
		strLower.push_back(tolower(strInput[i]));
		i++;
	}
	return strLower;
}
