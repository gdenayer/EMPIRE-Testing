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
 * It is coupled with meshClientB.
 * A dummy client which reads gid mesh, sends constant data field to the server.
 * If time step is 3, then in loose coupling, the data field sent is (3.0, 3.0, 3.0),
 * or in iterative coupling, the data field sent is (3.01, 3.01, 3.01) if meanwhile it is
 * at the 2nd iterative loop.
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
    int numNodes = -1;
    int numElems = -1;
    double *nodeCoors;
    int *nodeIDs;
    int *numNodesPerElem;
    int *elemTable;
    int *elemIDs;
    GiDFileIO::readDotMsh(meshFile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable, elemIDs);
    EMPIRE_API_sendMesh("defaultMesh", numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable);

    string resFile = EMPIRE_API_getUserDefinedText("GiDResFile");
    string dataFieldName = EMPIRE_API_getUserDefinedText("dataFieldName");
    string analysisName = EMPIRE_API_getUserDefinedText("analysisName");
    string stepNumStr = EMPIRE_API_getUserDefinedText("stepNum");
    int stepNum = getInteger(stepNumStr);

    string dataFieldNameDisp = "\"" + dataFieldName + "_disp\"";
    string dataFieldNameRot = "\"" + dataFieldName + "_rot\"";

    /*double RADIUS = 5.0;
     //double RADIUS = 10.0 / M_PI;
     // define deformation of curve/beam in Q
     double *dofsSend = new double[6 * numNodes];
     for (int i = 0; i < numNodes; i++) {
     double angle = 0.0 + (double) (i) * (M_PI / numElems);
     dofsSend[i * 6 + 0] = RADIUS * sin(angle) - nodeCoors[i * 3 + 0];
     dofsSend[i * 6 + 1] = (RADIUS - RADIUS * cos(angle)) / sqrt(2.0)
     - nodeCoors[i * 3 + 1];
     dofsSend[i * 6 + 2] = -(RADIUS - RADIUS * cos(angle)) / sqrt(2.0);
     dofsSend[i * 6 + 3] = 0.0;
     dofsSend[i * 6 + 4] = angle / sqrt(2.0);
     dofsSend[i * 6 + 5] = angle / sqrt(2.0);
     // this rotation vector cannot be transformed to between non-parallel systems
     }*/

    double *disp = new double[numNodes * 3];
    double *rot = new double[numNodes * 3];
    double *dofsSend = new double[numNodes * 6];
    double *dofsRecv = new double[numNodes * 6];
    for (int i = 1; i <= numTimeSteps; i++) {
        do {
            /*for (int j = 0; j < numNodes; j++)
             dofsSend[j * 3 + 2] = i;*/
            GiDFileIO::readNodalDataFromDotRes(resFile, dataFieldNameDisp, analysisName, stepNum,
                    "vector", numNodes, nodeIDs, disp);
            GiDFileIO::readNodalDataFromDotRes(resFile, dataFieldNameRot, analysisName, stepNum,
                    "vector", numNodes, nodeIDs, rot);
            for (int j = 0; j < numNodes; j++) {
                dofsSend[j * 6 + 0] = disp[j * 3 + 0];
                dofsSend[j * 6 + 1] = disp[j * 3 + 1];
                dofsSend[j * 6 + 2] = disp[j * 3 + 2];
                dofsSend[j * 6 + 3] = rot[j * 3 + 0];
                dofsSend[j * 6 + 4] = rot[j * 3 + 1];
                dofsSend[j * 6 + 5] = rot[j * 3 + 2];
            }

            for (int j = 0; j < numNodes * 6; j++)
                dofsSend[j] *= (double) i / (double) numTimeSteps;
            cout << "sending ..." << endl;
            EMPIRE_API_sendDataField("defaultField", numNodes * 6, dofsSend);
            //EMPIRE_API_printDataField("-----send-----", numNodes * 3, dofsSend);

            cout << "receiving ..." << endl;
            EMPIRE_API_recvDataField("defaultField", numNodes * 6, dofsRecv);
            //EMPIRE_API_printDataField("-----receive-----", numNodes * 3, dofsRecv);

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
    delete[] elemTable;
    delete[] numNodesPerElem;
    delete[] elemIDs;
    delete[] rot;
    delete[] disp;
    delete[] dofsSend;
    delete[] dofsRecv;
}
void standardCoupling() {
    bool todoIterativeCoupling = doIterativeCoupling(EMPIRE_API_getUserDefinedText("couplingType"));
    string meshFile = EMPIRE_API_getUserDefinedText("GiDMeshFile");
    int numTimeSteps = getInteger(EMPIRE_API_getUserDefinedText("numTimeSteps"));

    int numNodes = -1;
    int numElems = -1;
    double *nodeCoors;
    int *nodeIDs;
    int *numNodesPerElem;
    int *elemTable;
    int *elemIDs;
    GiDFileIO::readDotMsh(meshFile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable, elemIDs);
    EMPIRE_API_sendMesh("defaultMesh", numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable);

    double *dofsSend = new double[numNodes * 3];
    for (int j = 0; j < numNodes * 3; j++)
        dofsSend[j] = 0.0;
    double *dofsRecv = new double[numNodes * 3];
    for (int j = 0; j < numNodes * 3; j++)
        dofsRecv[j] = 0.0;

    FieldCreator *fc = new FieldCreator(numNodes, numElems, numNodesPerElem[0], nodeCoors, nodeIDs,
            elemTable, "constant");
    for (int i = 1; i <= numTimeSteps; i++) {
        if (!todoIterativeCoupling) {
            /*for (int j = 0; j < numNodes; j++)
             dofsSend[j * 3 + 2] = i;*/
            fc->create(dofsSend);
            for (int j = 0; j < numNodes * 3; j++)
                dofsSend[j] *= i;
            cout << "sending ..." << endl;
            EMPIRE_API_sendDataField("defaultField", numNodes * 3, dofsSend);
            //EMPIRE_API_printDataField("-----send-----", numNodes * 3, dofsSend);

            cout << "receiving ..." << endl;
            EMPIRE_API_recvDataField("defaultField", numNodes * 3, dofsRecv);
            //EMPIRE_API_printDataField("-----receive-----", numNodes * 3, dofsRecv);

            dummyUpdataPDE();
            dummySolvePDE();
        } else {
            double solution = (double) i;
            for (int j = 0; j < numNodes; j++)
                dofsSend[j * 3 + 2] = solution;

            int count = 0;
            do {
                count++;
                for (int j = 0; j < numNodes; j++)
                    dofsSend[j * 3 + 2] = solution + pow(10, (double) (-count));
                cout << "sending ..." << endl;
                EMPIRE_API_sendDataField("defaultField", numNodes * 3, dofsSend);
                //EMPIRE_API_printDataField("-----send-----", numNodes * 3, dofsSend);
                cout << "receiving ..." << endl;
                EMPIRE_API_recvDataField("defaultField", numNodes * 3, dofsRecv);
                //EMPIRE_API_printDataField("-----receive-----", numNodes * 3, dofsRecv);
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

