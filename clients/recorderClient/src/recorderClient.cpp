#include "EMPIRE_API.h"
#include "EMPEROR_Enum.h"
#include <assert.h>
#include "GiDFileIO.h"
#include "GiDIGAFileIO.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <cstdlib>

/* enum EMPIRE_Mesh_type {
    EMPIRE_Mesh_FEMesh,
    EMPIRE_Mesh_IGAMesh,
    EMPIRE_Mesh_copyFEMesh,
    EMPIRE_Mesh_SectionMesh
}; */

using namespace std;

void convertToLowerCase(string &str);

bool doIterativeCoupling(string text) {
    convertToLowerCase(text);
    if (text == "iterativecoupling")
        return true;
    else if (text == "loosecoupling")
        return false;
    else
        assert(false);
    return false;
}

bool onNodes(string text) {
    convertToLowerCase(text);
    if (text == "onnodes")
        return true;
    else if (text == "ongausspoints")
        return false;
    else
        assert(false);
    return false;
}

int getNumTimeSteps(string text) {
    assert(text != "");
    stringstream ss(text);
    int numTimeSteps;
    ss >> numTimeSteps;
    return numTimeSteps;
}

int getReadInterval(string text) {
    assert(text != "");
    stringstream ss(text);
    int readInterval;
    ss >> readInterval;
    return readInterval;
}

int getReadStart(string text) {
    assert(text != "");
    stringstream ss(text);
    int readStart;
    ss >> readStart;
    return readStart;
}

/*
 * RecorderClient reads GiD .msh and .res file, send result at each timestep to the server
 */
int main(int argc, char **argv) {
    // EMPIRE_API_Connect("recorderClientInput.xml");
    if (argc != 2) {
        std::cout << "Provide a valid input file" << std::endl;
        exit(-1);
    } else
        EMPIRE_API_Connect(argv[1]);

    int numTimeSteps = getNumTimeSteps(EMPIRE_API_getUserDefinedText("numTimeSteps"));
    int readInterval = getReadInterval(EMPIRE_API_getUserDefinedText("readInterval"));
    int readStart = getReadStart(EMPIRE_API_getUserDefinedText("readStart"));
    string meshType = EMPIRE_API_getUserDefinedText("meshType");
    string meshfile = EMPIRE_API_getUserDefinedText("GiDMeshFile");
    string resfileName = EMPIRE_API_getUserDefinedText("GiDResultFile");
    string resultName = EMPIRE_API_getUserDefinedText("GiDResultFile_resultName");
    string analysisName = EMPIRE_API_getUserDefinedText("GiDResultFile_analysisName");
    bool onNodes = EMPIRE_API_getUserDefinedText("GiDResultFile_dataFieldLocation");

    // Check the input mesh
    bool isFEMesh;
    if (meshType.compare("FEMesh")) {
        isFEMesh = false;
    } else if (meshType.compare("IGAMesh")) {
        isFEMesh = true;
    } else
        assert(false);

    int numNodes = -1;
    int numElems = -1;
    double *nodeCoors;
    int *nodeIDs;
    int *numNodesPerElem;
    int *elemTable;
    int *elemIDs;
    GiDFileIO::readDotMsh(meshfile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable, elemIDs);
    EMPIRE_API_sendMesh("defaultMesh", numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable);
    /* If send section mesh
    double rotation[9];
    for (int i = 0; i < 9; i++)
        rotation[i] = 0.0;
    rotation[0] = 1.0;
    rotation[4] = 1.0;
    rotation[8] = 1.0;
    double translation[3];
    for (int i = 0; i < 3; i++)
        translation[i] = 0.0;
    EMPIRE_API_sendSectionMesh("defaultMesh", numNodes, numElems, nodeCoors, nodeIDs,
            numNodesPerElem, elemTable, 20, 10, 10, 10, rotation, translation);*/

    bool todoIterativeCoupling = doIterativeCoupling(EMPIRE_API_getUserDefinedText("couplingType"));

    double *nodalData = new double[numNodes * 3];
    double *elementalData = new double[numElems * 3];

    ifstream dotResFile(resfileName.c_str());

    for (int i = 0; i < numTimeSteps; i++) {
        std::cout << "Read time step " << readStart + i * readInterval << " in .res" << std::endl;
        do {
            if (onNodes) {
                GiDFileIO::readNodalDataFromDotResFast(dotResFile, resfileName, resultName, analysisName,
                        readStart + i * readInterval, "vector", numNodes, nodeIDs, nodalData);
                EMPIRE_API_sendDataField("defaultField", numNodes * 3, nodalData);
            } else {
                GiDFileIO::readElementalDataFromDotResFast(dotResFile, resfileName, resultName, analysisName,
                        readStart + i * readInterval, "vector", numElems, elemIDs, numNodesPerElem,
                        elementalData);
                EMPIRE_API_sendDataField("defaultField", numElems * 3, elementalData);
            }

            // do not receive any data
            if (!todoIterativeCoupling)
                break;
        } while (EMPIRE_API_recvConvergenceSignal() == 0);
    }

    if (meshType.compare("FEMMesh")) {
        delete[] nodalData;
        delete[] elementalData;
        delete[] nodeCoors;
        delete[] nodeIDs;
        delete[] numNodesPerElem;
        delete[] elemTable;
        delete[] elemIDs;
    }

    EMPIRE_API_Disconnect();
    return (0);
}

/***********************************************************************************************
 * \brief Convert the string to a lower case string
 * \param[in] str the string
 * \author Tianyang Wang
 ***********/
void convertToLowerCase(string &str) {
    for (unsigned i = 0; i < str.size(); i++) {
        str[i] = tolower(str[i]);
    }
}
