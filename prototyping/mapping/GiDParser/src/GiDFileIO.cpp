#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <assert.h>

#include "GiDFileIO.h"
#include "GiDBasicData.h"
#include "GiDResultInventory.h"
#include "GiDMeshInventory.h"

using namespace std;

// a static function means the function can only be called in this file!!!
static void GiDFileIO_extractGiDMesh(GiDMeshInventory *MI, int &numNodes, int &numElems,
        int &nodesPerElem, double *&nodeCoors, int *&nodeNumbers, int *&elemTable);

void GiDFileIO::readDotMsh(string fileName, int &numNodes, int &numElems, int &nodesPerElem,
        double *&nodeCoors, int *&nodeNumbers, int *&elemTable) {
    GiDMeshInventory *MI = new GiDMeshInventory(fileName);
    GiDFileIO_extractGiDMesh(MI, numNodes, numElems, nodesPerElem, nodeCoors, nodeNumbers,
            elemTable);
    delete MI;
}

void GiDFileIO::writeDotMsh(std::string fileName, int numNodes, int numElems, int nodesPerElem,
        double *nodeCoors, int *nodeNumbers, int *elemTable) {
    vector<GiDNode*> *nodes = new vector<GiDNode*>;
    for (int i = 0; i < numNodes; i++) {
        double *coor = new double[3];
        for (int j = 0; j < 3; j++)
            coor[j] = nodeCoors[i * 3 + j];
        GiDNode *node = new GiDNode(nodeNumbers[i], coor);
        nodes->push_back(node);
    }
    vector<GiDElement*> *elements = new vector<GiDElement*>;
    for (int i = 0; i < numElems; i++) {
        int *elem = new int[nodesPerElem];
        for (int j = 0; j < nodesPerElem; j++)
            elem[j] = elemTable[i * nodesPerElem + j];
        GiDElement *element = new GiDElement(i + 1, elem);
        elements->push_back(element);
    }
    string meshName = fileName;
    meshName.erase(meshName.size() - 4, 4);
    vector<GiDMesh*> *meshes = new vector<GiDMesh*>;
    string elementType;
    if (nodesPerElem == 3)
        elementType = "Triangle";
    else if (nodesPerElem == 4)
        elementType = "Quadrilateral";
    else
        assert(false);

    meshes->push_back(new GiDMesh(meshName, elementType, 3, nodesPerElem, nodes, elements));
    GiDMeshInventory *MI = new GiDMeshInventory(meshes);
    ofstream meshFile(fileName.c_str());
    MI->writeDotMsh(meshFile);
}

void GiDFileIO::writeDotRes(string fileName, int numNodes, int *nodeNumbers, string analysisName,
        int numResults, string *names, double **results) {

    vector<GiDResult*> *vec2 = new vector<GiDResult*>();
    GiDResultInventory *RI = new GiDResultInventory(vec2);

    string type = "Vector";
    int stepNum = 1;
    for (int k = 0; k < numResults; k++) {
        vector<GiDValue*> *vec = new vector<GiDValue*>();
        const int DOF_DIM = 3; // only allow vector!!!
        for (int i = 0; i < numNodes; i++) {
            double* dofs = new double[DOF_DIM];
            for (int j = 0; j < DOF_DIM; j++)
                dofs[j] = results[k][i * DOF_DIM + j];
            int id = nodeNumbers[i];
            vec->push_back(new GiDValue(id, dofs));
        }

        GiDResult *res = new GiDResult(names[k], analysisName, stepNum, type, vec);
        RI->getResults()->push_back(res);
    }
    ofstream outFile(fileName.c_str());
    cout << "Write DOFs to " << fileName << endl;
    RI->writeDotRes(outFile);
    delete RI;
}

GiDResultInventory *GiDFileIO::initDotRes() {
    vector<GiDResult*> *results = new vector<GiDResult*>();
    return new GiDResultInventory(results);
}

void GiDFileIO::appendResultToDotRes(GiDResultInventory *ri, std::string resultName,
        std::string analysisName, int stepNum, std::string type, int numNodes, int *nodeNumbers,
        double *data) {

    assert(ri != NULL);

    int dimension = 0;
    if (type == "Vector")
        dimension = 3;
    else if (type == "Scalar")
        dimension = 1;
    else
        assert(false);

    vector<GiDValue*> *values = new vector<GiDValue*>();
    for (int i = 0; i < numNodes; i++) {
        double* nodalDof = new double[dimension];
        for (int j = 0; j < dimension; j++)
            nodalDof[j] = data[i * dimension + j];
        int id = nodeNumbers[i];
        values->push_back(new GiDValue(id, nodalDof));
    }

    GiDResult *res = new GiDResult(resultName, analysisName, stepNum, type, values);
    ri->getResults()->push_back(res);
}

void GiDFileIO::writeDotRes(GiDResultInventory *ri, std::string fileName) {
    assert(ri != NULL);
    ofstream outFile(fileName.c_str());
    cout << "Write DOFs to " << fileName << endl;
    ri->writeDotRes(outFile);
    delete ri;
}

void GiDFileIO::writeGnuplotError(string testCaseName, int numNodes, double *nodeCoors,
        double *errorVec) {
    const string GNUPLOT = "_error.gnuplot";
    const int XYZ = 3;
    string fileName = testCaseName;
    fileName.append(GNUPLOT);
    cout << "write DOFs to " << fileName << endl;
    ofstream outFile(fileName.c_str());
    outFile << "# output by GiDFileIO" << endl;
    for (int i = 0; i < numNodes; i++) {
        double errorNorm = 0.0;
        errorNorm += errorVec[i * XYZ + 0] * errorVec[i * XYZ + 0];
        errorNorm += errorVec[i * XYZ + 1] * errorVec[i * XYZ + 1];
        errorNorm += errorVec[i * XYZ + 2] * errorVec[i * XYZ + 2];
        errorNorm = sqrt(errorNorm);
        outFile << nodeCoors[i * XYZ] << '\t' << nodeCoors[i * XYZ + 1] << '\t'
                << nodeCoors[i * XYZ + 2] << '\t' << errorNorm << endl;
    }
}

void GiDFileIO_extractGiDMesh(GiDMeshInventory *MI, int &numNodes, int &numElems, int &nodesPerElem,
        double *&nodeCoors, int *&nodeNumbers, int *&elemTable) {
    const int XYZ = 3;
    numNodes = MI->getMeshes()->at(0)->getNodes()->size();
    numElems = MI->getMeshes()->at(0)->getElements()->size();
    nodesPerElem = MI->getMeshes()->at(0)->getNodesPerElem();

    nodeCoors = new double[numNodes * XYZ];
    nodeNumbers = new int[numNodes];
    elemTable = new int[numElems * nodesPerElem];

    for (int i = 0; i < numNodes; i++)
        for (int j = 0; j < XYZ; j++)
            nodeCoors[i * XYZ + j] = MI->getMeshes()->at(0)->getNodes()->at(i)->getCoordinate()[j];

    for (int i = 0; i < numNodes; i++)
        nodeNumbers[i] = MI->getMeshes()->at(0)->getNodes()->at(i)->getId();

    for (int i = 0; i < numElems; i++)
        for (int j = 0; j < nodesPerElem; j++)
            elemTable[i * nodesPerElem + j] =
                    MI->getMeshes()->at(0)->getElements()->at(i)->getNodeIds()[j];
}

