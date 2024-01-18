#include "EMPIRE_API.h"
#include <assert.h>
#include "GiDFileIO.h"
#include "FluxCreator.h"
#include "FieldCreator.h"
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

void dummyUpdataPDE() {
    // t = t + delta_t;
    // set up the PDE at this time step
}
void dummySolvePDE() {
    // solve the PDE
}

enum Scheme {
    Gauss_Seidel, Jacobi
} scheme = Jacobi;

bool doIterativeCoupling(string text) {
    if (text == "iterativeCoupling")
        return true;
    else if (text == "looseCoupling")
        return false;
    else
        assert(false);
    return false;
}

int getInteger(string text) {
    stringstream ss(text);
    int i;
    ss >> i;
    return i;
}

void curveSurfaceMapping();
void standardCoupling();
/*
 * It is coupled with meshClientA.
 * A dummy client which reads gid mesh, sends consistent nodal forces to the server.
 * The consistent nodal force is computed from a constant pressure field over the surface.
 */
int main(int argc, char **argv) {
    if (argc != 2) {
        cerr << "Error: Please provide a valid input file." << endl;
        exit (EXIT_FAILURE);
    }
    EMPIRE_API_Connect(argv[1]);

    string isCurveSurfaceMapping(EMPIRE_API_getUserDefinedText("curveSurfaceMapping"));
    if (isCurveSurfaceMapping == "true") {
        curveSurfaceMapping();
    } else {
        standardCoupling();
    }

    EMPIRE_API_Disconnect();
    return (0);
}

void curveSurfaceMapping() {
    bool todoIterativeCoupling = doIterativeCoupling(EMPIRE_API_getUserDefinedText("couplingType"));
    int numTimeSteps = getInteger(EMPIRE_API_getUserDefinedText("numTimeSteps"));

    string meshFile = EMPIRE_API_getUserDefinedText("GiDMeshFile");
    // surface mesh in Q
    int numNodes = -1;
    int numElems = -1;
    double *nodeCoors;
    int *nodeIDs;
    int *numNodesPerElem;
    int *elemTable;
    int *elemIDs;
    GiDFileIO::readDotMsh(meshFile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable, elemIDs);

    // section information
    string numSectionsString = EMPIRE_API_getUserDefinedText("numSections");
    string numRootSectionNodesString = EMPIRE_API_getUserDefinedText("numRootSectionNodes");
    string numNormalSectionNodesString = EMPIRE_API_getUserDefinedText("numNormalSectionNodes");
    string numTipSectionNodesString = EMPIRE_API_getUserDefinedText("numTipSectionNodes");
    int numSections = getInteger(numSectionsString);
    int numRootSectionNodes = getInteger(numRootSectionNodesString);
    int numNormalSectionNodes = getInteger(numNormalSectionNodesString);
    int numTipSectionNodes = getInteger(numTipSectionNodesString);

    // compute km
    double rotation[9];
    for (int i = 0; i < 9; i++)
        rotation[i] = 0.0;
    rotation[0] = 1.0;
    rotation[4] = 1.0;
    rotation[8] = 1.0;
    double translation[3];
    for (int i = 0; i < 3; i++)
        translation[i] = 0.0;
    // disp of surface
    double *dofsRecv = new double[numNodes * 3];
    double *dofsSend = new double[numNodes * 3];

    EMPIRE_API_sendSectionMesh("defaultMesh", numNodes, numElems, nodeCoors, nodeIDs,
            numNodesPerElem, elemTable, numSections, numRootSectionNodes, numNormalSectionNodes,
            numTipSectionNodes, rotation, translation);

    for (int i = 1; i <= numTimeSteps; i++) {
        do {
            cout << "receiving ..." << endl;
            EMPIRE_API_recvDataField("defaultField", numNodes * 3, dofsRecv);
            //EMPIRE_API_printDataField("-----receive-----", numNodes * 3, dofsRecv);
            /*for (int j = 0; j < numNodes; j++) {
                dofsSend[j * 3 + 0] = 0.0;
                dofsSend[j * 3 + 1] = 50.0;
                dofsSend[j * 3 + 2] = 0.0;
            }*/
            for (int j = 0; j < numNodes; j++) {
                dofsSend[j * 3 + 0] = 0.0;
                dofsSend[j * 3 + 1] = 0.0;
                dofsSend[j * 3 + 2] = 0.0;
            }
            //dofsSend[12 * 3 + 1] = 10000.0;
            cout << "sending ..." << endl;
            EMPIRE_API_sendDataField("defaultField", numNodes * 3, dofsSend);
            //EMPIRE_API_printDataField("-----send-----", numNodes * 3, dofsSend);
            if (!todoIterativeCoupling) {
                break;
            }
        } while (EMPIRE_API_recvConvergenceSignal() == 0);
        /*std::cout << std::endl;
         std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
         << std::endl;*/
    }
    delete[] nodeCoors;
    delete[] nodeIDs;
    delete[] numNodesPerElem;
    delete[] elemTable;
    delete[] elemIDs;
    delete[] dofsRecv;
    delete[] dofsSend;
}

void standardCoupling() {
    bool todoIterativeCoupling = doIterativeCoupling(EMPIRE_API_getUserDefinedText("couplingType"));
    string meshFile = EMPIRE_API_getUserDefinedText("GiDMeshFile");
    int numTimeSteps = getInteger(EMPIRE_API_getUserDefinedText("numTimeSteps"));

    int numNodes;
    int numElems;
    double *nodeCoors;
    int *nodeIDs;
    int *numNodesPerElem;
    int *elemTable;
    int *elemIDs;

    GiDFileIO::readDotMsh(meshFile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable, elemIDs);

    EMPIRE_API_sendMesh("defaultMesh", numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable);

    FluxCreator *fc = new FluxCreator(numNodes, numElems, numNodesPerElem[0], nodeCoors, nodeIDs,
            elemTable, "constant", true);

    bool sendTractionElem = false;
    string sendDataField = EMPIRE_API_getUserDefinedText("sendDataField");
    if (sendDataField == "traction")
        sendTractionElem = true;
    else if (sendDataField == "force")
        sendTractionElem = false;
    else
        assert(false);

    double *dofsSend;
    if (!sendTractionElem) {
        dofsSend = new double[numNodes * 3];
        for (int j = 0; j < numNodes * 3; j++)
            dofsSend[j] = 0.0;
    } else {
        dofsSend = new double[numElems * 3];
        for (int j = 0; j < numElems * 3; j++)
            dofsSend[j] = 0.0;
    }
    double *dofsRecv = new double[numNodes * 3];
    for (int j = 0; j < numNodes * 3; j++)
        dofsRecv[j] = 0.0;

    for (int i = 1; i <= numTimeSteps; i++) {
        if (!todoIterativeCoupling) {
            cout << "receiving ..." << endl;
            EMPIRE_API_recvDataField("defaultField", numNodes * 3, dofsRecv);
            //EMPIRE_API_printDataField("-----receive-----", numNodes * 3, dofsRecv);

            if (scheme == Gauss_Seidel) {
                dummyUpdataPDE();
                dummySolvePDE();
            }

            if (!sendTractionElem) {
                fc->create(dofsSend);
                for (int j = 0; j < numNodes * 3; j++)
                    dofsSend[j] *= i;
                cout << "sending ..." << endl;
                EMPIRE_API_sendDataField("defaultField", numNodes * 3, dofsSend);
                //EMPIRE_API_printDataField("-----send-----", numNodes * 3, dofsSend);
            } else {
                for (int j = 0; j < numElems; j++)
                    dofsSend[j * 3 + 2] = i;
                cout << "sending ..." << endl;
                EMPIRE_API_sendDataField("defaultField", numElems * 3, dofsSend);
                //EMPIRE_API_printDataField("-----send-----", numElems * 3, dofsSend);
            }

            if (scheme == Jacobi) {
                dummyUpdataPDE();
                dummySolvePDE();
            }
        } else {
            do {
                cout << "receiving ..." << endl;
                EMPIRE_API_recvDataField("defaultField", numNodes * 3, dofsRecv);
                //EMPIRE_API_printDataField("-----receive-----", numNodes * 3, dofsRecv);

                if (!sendTractionElem) {
                    fc->create(dofsSend);
                    for (int j = 0; j < numNodes; j++)
                        dofsSend[j * 3 + 2] *= i;
                    cout << "sending ..." << endl;
                    EMPIRE_API_sendDataField("defaultField", numNodes * 3, dofsSend);
                    //EMPIRE_API_printDataField("-----send-----", numNodes * 3, dofsSend);
                } else {
                    for (int j = 0; j < numElems; j++)
                        dofsSend[j * 3 + 2] = i;
                    cout << "sending ..." << endl;
                    EMPIRE_API_sendDataField("defaultField", numElems * 3, dofsSend);
                    //EMPIRE_API_printDataField("-----send-----", numElems * 3, dofsSend);
                }
            } while (EMPIRE_API_recvConvergenceSignal() == 0);
            /*std::cout << std::endl;
             std::cout << "Time step: " << i << ", no. of inner loop steps: " << count
             << std::endl;*/
        }
    }

    delete[] nodeCoors;
    delete[] nodeIDs;
    delete[] numNodesPerElem;
    delete[] elemTable;
    delete[] elemIDs;
    delete[] dofsRecv;
    delete[] dofsSend;
}

