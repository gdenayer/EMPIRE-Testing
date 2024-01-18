#include "EMPIRE_API.h"
#include <assert.h>
#include "GiDFileIO.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>

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

int getNumTimeSteps(string text) {
    assert(text != "");
    stringstream ss(text);
    int numTimeSteps;
    ss >> numTimeSteps;
    return numTimeSteps;
}

double getTimeInterval(string text) {
    assert(text != "");
    stringstream ss(text);
    double timeInterval;
    ss >> timeInterval;
    return timeInterval;
}

double getOmega(string text) {
    assert(text != "");
    stringstream ss(text);
    double omega;
    ss >> omega;
    return omega;
}

void computeGravityForces(double numNodes, double *XX, double *MASSPerLength, double *forces) {
    for (int i = 0; i < numNodes * 6; i++) {
        forces[i] = 0.0;
    }
    for (int i = 0; i < numNodes - 1; i++) {
        double f = (MASSPerLength[i] + MASSPerLength[i + 1]) / 2.0 * fabs(XX[i + 1] - XX[i]) * (-9.8);
        forces[i * 6 + 2] += f / 2.0;
        forces[(i + 1) * 6 + 2] += f / 2.0;
    }
}
void computeCentrifugalForces(double omega, double numNodes, double *XX, double *MASSPerLength,
        double *forces) {
    for (int i = 0; i < numNodes * 6; i++) {
        forces[i] = 0.0;
    }

    for (int i = 0; i < numNodes - 1; i++) {
        double f = (MASSPerLength[i] + MASSPerLength[i + 1]) / 2.0 * fabs(XX[i + 1] - XX[i]) * omega
                * omega * (XX[i + 1] + XX[i]) / 2.0; // sign of XX decides the direction of the force
        forces[i * 6 + 0] += f / 2.0;
        forces[(i + 1) * 6 + 0] += f / 2.0;
    }

}
void rotateForces(double phi, int numNodes, double *forces) {
    phi = -phi;
    for (int i = 0; i < numNodes; i++) {
        double f_x = cos(phi) * forces[i * 6 + 0] - sin(phi) * forces[i * 6 + 2];
        double f_z = sin(phi) * forces[i * 6 + 0] + cos(phi) * forces[i * 6 + 2];
        forces[i * 6 + 0] = f_x;
        forces[i * 6 + 2] = f_z;
        double m_x = cos(phi) * forces[i * 6 + 3] - sin(phi) * forces[i * 6 + 5];
        double m_z = sin(phi) * forces[i * 6 + 3] + cos(phi) * forces[i * 6 + 5];
        forces[i * 6 + 3] = m_x;
        forces[i * 6 + 5] = m_z;
    }
}

/*
 * NRELBladeBeam reads GiD .msh and .res file, send result at each timestep to the server
 */
int main(int argc, char **argv) {
    EMPIRE_API_Connect("NRELBladeBeam.xml");

    int numTimeSteps = getNumTimeSteps(EMPIRE_API_getUserDefinedText("numTimeSteps"));
    double omega = getOmega(EMPIRE_API_getUserDefinedText("omega"));
    double timeInterval = getTimeInterval(EMPIRE_API_getUserDefinedText("timeInterval"));
    bool todoIterativeCoupling = doIterativeCoupling(EMPIRE_API_getUserDefinedText("couplingType"));

    int numberOfNodes;
    double *MASSPerLength;
    double *XX_plus;
    double *XX_minus;
    { // blade in X+
        int numNodes = -1;
        int numElems = -1;
        double *nodeCoors;
        int *nodeIDs;
        int *numNodesPerElem;
        int *elemTable;
        int *elemIDs;
        string meshFilePlus = EMPIRE_API_getUserDefinedText("GiDMeshFilePlus");
        GiDFileIO::readDotMsh(meshFilePlus, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
                elemTable, elemIDs);
        EMPIRE_API_sendMesh("", numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable);

        numberOfNodes = numNodes;
        XX_plus = new double[numNodes];
        for (int i = 0; i < numNodes; i++) {
            XX_plus[i] = nodeCoors[i * 3 + 0];
        }

        string resFilePlus = EMPIRE_API_getUserDefinedText("GiDResFilePlus");
        MASSPerLength = new double[numNodes];
        GiDFileIO::readNodalDataFromDotRes(resFilePlus, "MASSPerLength", "EMPIRE_cosimulation", 1, "scalar",
                numNodes, nodeIDs, MASSPerLength);

        delete[] nodeCoors;
        delete[] nodeIDs;
        delete[] numNodesPerElem;
        delete[] elemTable;
        delete[] elemIDs;
    }
    { // blade in X-
        int numNodes = -1;
        int numElems = -1;
        double *nodeCoors;
        int *nodeIDs;
        int *numNodesPerElem;
        int *elemTable;
        int *elemIDs;
        string meshFile = EMPIRE_API_getUserDefinedText("GiDMeshFileMinus");
        GiDFileIO::readDotMsh(meshFile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
                elemTable, elemIDs);
        EMPIRE_API_sendMesh("", numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable);

        XX_minus = new double[numNodes];
        for (int i = 0; i < numNodes; i++) {
            XX_minus[i] = nodeCoors[i * 3 + 0];
        }

        delete[] nodeCoors;
        delete[] nodeIDs;
        delete[] numNodesPerElem;
        delete[] elemTable;
        delete[] elemIDs;
    }

    double centrifugalForce_plus[numberOfNodes * 6];
    double centrifugalForce_minus[numberOfNodes * 6];
    double gravityForce_plus[numberOfNodes * 6];
    double gravityForce_minus[numberOfNodes * 6];
    double forces_plus[numberOfNodes * 6];
    double forces_minus[numberOfNodes * 6];
    computeCentrifugalForces(omega, numberOfNodes, XX_plus, MASSPerLength, centrifugalForce_plus);
    computeCentrifugalForces(omega, numberOfNodes, XX_minus, MASSPerLength, centrifugalForce_minus);
    computeGravityForces(numberOfNodes, XX_plus, MASSPerLength, gravityForce_plus);
    computeGravityForces(numberOfNodes, XX_minus, MASSPerLength, gravityForce_minus);

    for (int i = 0; i < numTimeSteps; i++) {
        std::cout << "Do time step " << i + 1 << " ..." << std::endl;
        double phi = omega * timeInterval * (double) (i + 1);
        do {
            EMPIRE_API_recvDataField("defaultField", numberOfNodes * 6, forces_plus);
            EMPIRE_API_recvDataField("defaultField", numberOfNodes * 6, forces_minus);
            for (int j = 0; j < numberOfNodes * 6; j++) {
                forces_plus[j] += gravityForce_plus[j];
                forces_minus[j] += gravityForce_minus[j];
            }
            for (int j = 0; j < numberOfNodes; j++) {
                forces_plus[j*6+3] = 0.0;
                forces_plus[j*6+4] = 0.0;
                forces_plus[j*6+5] = 0.0;
                forces_minus[j*6+3] = 0.0;
                forces_minus[j*6+4] = 0.0;
                forces_minus[j*6+5] = 0.0;
            }
            rotateForces(phi, numberOfNodes, forces_plus);
            rotateForces(phi, numberOfNodes, forces_minus);
            for (int j = 0; j < numberOfNodes * 6; j++) {
                forces_plus[j] += centrifugalForce_plus[j];
                forces_minus[j] += centrifugalForce_minus[j];
            }
            EMPIRE_API_sendDataField("defaultField", numberOfNodes * 6, forces_plus);
            EMPIRE_API_sendDataField("defaultField", numberOfNodes * 6, forces_minus);

            // do not receive any data
            if (!todoIterativeCoupling)
                break;
        } while (EMPIRE_API_recvConvergenceSignal() == 0);
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
