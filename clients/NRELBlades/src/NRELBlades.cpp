#ifdef USE_INTEL_MKL
#include <mkl.h>
#include <mkl_lapacke.h>
#include <mkl_spblas.h>
#endif

#ifndef USE_INTEL_MKL
#include <Dense>
#endif

#include "EMPIRE_API.h"
#include <assert.h>
#include "GiDFileIO.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include "FEMMath.h"
#include "Message.h"
#include "ConstantsAndVariables.h"
#include <assert.h>
#include <typeinfo>
#include <cmath>
#include "AuxiliaryParameters.h"

using namespace std;
using namespace EMPIRE::MathLibrary;

/********//**
 * \brief Class GaussQuadrature base class for general quadrature
 * \author Andreas Apostolatos
 ***********/
class GaussQuadrature {

private:
    // The number of Gauss Points
    int numGaussPoints;

    // Number of dimensions
    int numDimensions;

    // Array containing the Gauss Point locations in the quadrature space
    const double *gaussPoints;

    // Array containing the quadrature weights
    const double *weights;

    /// Constructor, destructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    GaussQuadrature(int _numGaussPoints, int _numDimensions) :
            numGaussPoints(_numGaussPoints), numDimensions(_numDimensions) {
    }

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~GaussQuadrature() {
    }

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Returns the number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    double getNumGaussPoints() {
        return numGaussPoints;
    }

    /***********************************************************************************************
     * \brief Returns the coordinates of the Gauss point
     * \param[in] _index Index to the Gauss point
     * \author Andreas Apostolatos
     ***********/
    const double* getGaussPoint(int _index) {
        return &gaussPoints[_index * numDimensions];
    }

    /***********************************************************************************************
     * \brief Returns the weight of the Gauss point
     * \param[in] _index Index to the Gauss weight
     * \author Andreas Apostolatos
     ***********/
    double getGaussWeight(int _index) {
        return weights[_index];
    }

    /***********************************************************************************************
     * \brief Sets the coordinates of the Gauss points
     * \param[in] _values Pointer to the array with the values
     * \author Andreas Apostolatos
     ***********/
    void setGaussPoints(const double* _values) {
        gaussPoints = _values;
    }

    /***********************************************************************************************
     * \brief Sets the Gauss weights
     * \param[in] _values Pointer to the array with the values
     * \author Andreas Apostolatos
     ***********/
    void setGaussWeights(const double* _values) {
        weights = _values;
    }
};

/********//**
 * \brief Class GaussQuadratureOnTriangle performs Gauss quadrature on triangle using a symmetric rule
 * \author Andreas Apostolatos
 ***********/
class GaussQuadratureTriangle: public GaussQuadrature {
public:
    /***********************************************************************************************
     * \brief Constructor
     * param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    GaussQuadratureTriangle(int _numGaussPoints);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~GaussQuadratureTriangle() {
    }
};

/********//**
 * \brief Class GaussQuadratureOnTriangle performs Gauss quadrature on canonical triangle using the degenerated quadrilateral
 * \author Andreas Apostolatos
 ***********/
class GaussQuadratureTriangleUsingDegeneratedQuadrilateral: public GaussQuadrature {
public:
    /***********************************************************************************************
     * \brief Constructor
     * param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    GaussQuadratureTriangleUsingDegeneratedQuadrilateral(int _numGaussPoints);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~GaussQuadratureTriangleUsingDegeneratedQuadrilateral() {
    }
};

/**********
 * \brief Class GaussQuadratureOnQuad performs Gauss quadrature on quad
 * \author Andreas Apostolatos, Altug Emiroglu
 ***********/
class GaussQuadratureBiunitInterval: public GaussQuadrature {
public:
    /***********************************************************************************************
     * \brief Constructor
     * param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    GaussQuadratureBiunitInterval(int _numGaussPoints);

    /***********************************************************************************************
     * \brief DestructorEmperor_INCLUDES
     * \author Andreas Apostolatos
     ***********/
    virtual ~GaussQuadratureBiunitInterval() {
    }
};

/**********
 * \brief Class GaussQuadratureOnQuad performs Gauss quadrature on quad
 * \author Andreas Apostolatos
 ***********/
class GaussQuadratureBiunitQuadrilateral: public GaussQuadrature {
public:
    /***********************************************************************************************
     * \brief Constructor
     * param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    GaussQuadratureBiunitQuadrilateral(int _numGaussPoints);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~GaussQuadratureBiunitQuadrilateral() {
    }
};

GaussQuadratureTriangle::GaussQuadratureTriangle(int _numGaussPoints) :
        GaussQuadrature(_numGaussPoints, 2) {
    switch (_numGaussPoints) {
    case 1:
        setGaussPoints(IGACanonicalTriangleGaussPoints1);
        setGaussWeights(IGACanonicalTriangleWeights1);
        break;
    case 3:
        setGaussPoints(IGACanonicalTriangleGaussPoints3);
        setGaussWeights(IGACanonicalTriangleWeights3);
        break;
    case 4:
        setGaussPoints(IGACanonicalTriangleGaussPoints4);
        setGaussWeights(IGACanonicalTriangleWeights4);
        break;
    case 6:
        setGaussPoints(IGACanonicalTriangleGaussPoints6);
        setGaussWeights(IGACanonicalTriangleWeights6);
        break;
    case 7:
        setGaussPoints(IGACanonicalTriangleGaussPoints7);
        setGaussWeights(IGACanonicalTriangleWeights7);
        break;
    case 12:
        setGaussPoints(IGACanonicalTriangleGaussPoints12);
        setGaussWeights(IGACanonicalTriangleWeights12);
        break;
    case 13:
        setGaussPoints(IGACanonicalTriangleGaussPoints13);
        setGaussWeights(IGACanonicalTriangleWeights13);
        break;
    case 16:
        setGaussPoints(IGACanonicalTriangleGaussPoints16);
        setGaussWeights(IGACanonicalTriangleWeights16);
        break;
    default:
        assert(false);
    }
}

GaussQuadratureTriangleUsingDegeneratedQuadrilateral::GaussQuadratureTriangleUsingDegeneratedQuadrilateral(int _numGaussPoints) :
        GaussQuadrature(_numGaussPoints, 2) {
    switch (_numGaussPoints) {
    case 1:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints1);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights1);
        break;
    case 4:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints4);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights4);
        break;
    case 9:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints9);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights9);
        break;
    case 16:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints16);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights16);
        break;
    case 25:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints25);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights25);
        break;
    case 36:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints36);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights36);
        break;
    case 49:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints49);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights49);
        break;
    case 64:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints64);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights64);
        break;
    case 81:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints81);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights81);
        break;
    case 100:
        setGaussPoints(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussPoints100);
        setGaussWeights(IGACanonicalTriangleUsingDegeneratedQuadrilateralGaussWeights100);
        break;
    default:
        assert(false);
    }
}

GaussQuadratureBiunitInterval::GaussQuadratureBiunitInterval(int _numGaussPoints) :
    GaussQuadrature(_numGaussPoints, 1) {
    switch (_numGaussPoints) {
    case 1:
        setGaussPoints(IGABiunitIntervalGaussPoints1);
        setGaussWeights(IGABiunitIntervalGaussWeights1);
        break;
    case 2:
        setGaussPoints(IGABiunitIntervalGaussPoints2);
        setGaussWeights(IGABiunitIntervalGaussWeights2);
        break;
    case 3:
        setGaussPoints(IGABiunitIntervalGaussPoints3);
        setGaussWeights(IGABiunitIntervalGaussWeights3);
        break;
    case 4:
        setGaussPoints(IGABiunitIntervalGaussPoints4);
        setGaussWeights(IGABiunitIntervalGaussWeights4);
        break;
    case 5:
        setGaussPoints(IGABiunitIntervalGaussPoints5);
        setGaussWeights(IGABiunitIntervalGaussWeights5);
        break;
    case 6:
        setGaussPoints(IGABiunitIntervalGaussPoints6);
        setGaussWeights(IGABiunitIntervalGaussWeights6);
        break;
    case 7:
        setGaussPoints(IGABiunitIntervalGaussPoints7);
        setGaussWeights(IGABiunitIntervalGaussWeights7);
        break;
    case 8:
        setGaussPoints(IGABiunitIntervalGaussPoints8);
        setGaussWeights(IGABiunitIntervalGaussWeights8);
        break;
    case 9:
        setGaussPoints(IGABiunitIntervalGaussPoints9);
        setGaussWeights(IGABiunitIntervalGaussWeights9);
        break;
    case 10:
        setGaussPoints(IGABiunitIntervalGaussPoints10);
        setGaussWeights(IGABiunitIntervalGaussWeights10);
        break;
    case 11:
        setGaussPoints(IGABiunitIntervalGaussPoints11);
        setGaussWeights(IGABiunitIntervalGaussWeights11);
        break;
    case 12:
        setGaussPoints(IGABiunitIntervalGaussPoints12);
        setGaussWeights(IGABiunitIntervalGaussWeights12);
        break;
    case 13:
        setGaussPoints(IGABiunitIntervalGaussPoints13);
        setGaussWeights(IGABiunitIntervalGaussWeights13);
        break;
    case 14:
        setGaussPoints(IGABiunitIntervalGaussPoints14);
        setGaussWeights(IGABiunitIntervalGaussWeights14);
        break;
    case 15:
        setGaussPoints(IGABiunitIntervalGaussPoints15);
        setGaussWeights(IGABiunitIntervalGaussWeights15);
        break;
    case 16:
        setGaussPoints(IGABiunitIntervalGaussPoints16);
        setGaussWeights(IGABiunitIntervalGaussWeights16);
        break;
    case 17:
        setGaussPoints(IGABiunitIntervalGaussPoints17);
        setGaussWeights(IGABiunitIntervalGaussWeights17);
        break;
    case 18:
        setGaussPoints(IGABiunitIntervalGaussPoints18);
        setGaussWeights(IGABiunitIntervalGaussWeights18);
        break;
    case 19:
        setGaussPoints(IGABiunitIntervalGaussPoints19);
        setGaussWeights(IGABiunitIntervalGaussWeights19);
        break;
    case 20:
        setGaussPoints(IGABiunitIntervalGaussPoints20);
        setGaussWeights(IGABiunitIntervalGaussWeights20);
        break;
    case 21:
        setGaussPoints(IGABiunitIntervalGaussPoints21);
        setGaussWeights(IGABiunitIntervalGaussWeights21);
        break;
    case 22:
        setGaussPoints(IGABiunitIntervalGaussPoints22);
        setGaussWeights(IGABiunitIntervalGaussWeights22);
        break;
    case 23:
        setGaussPoints(IGABiunitIntervalGaussPoints23);
        setGaussWeights(IGABiunitIntervalGaussWeights23);
        break;
    case 24:
        setGaussPoints(IGABiunitIntervalGaussPoints24);
        setGaussWeights(IGABiunitIntervalGaussWeights24);
        break;
    case 25:
        setGaussPoints(IGABiunitIntervalGaussPoints25);
        setGaussWeights(IGABiunitIntervalGaussWeights25);
        break;
    case 26:
        setGaussPoints(IGABiunitIntervalGaussPoints26);
        setGaussWeights(IGABiunitIntervalGaussWeights26);
        break;
    case 27:
        setGaussPoints(IGABiunitIntervalGaussPoints27);
        setGaussWeights(IGABiunitIntervalGaussWeights27);
        break;
    case 28:
        setGaussPoints(IGABiunitIntervalGaussPoints28);
        setGaussWeights(IGABiunitIntervalGaussWeights28);
        break;
    case 29:
        setGaussPoints(IGABiunitIntervalGaussPoints29);
        setGaussWeights(IGABiunitIntervalGaussWeights29);
        break;
    case 30:
        setGaussPoints(IGABiunitIntervalGaussPoints30);
        setGaussWeights(IGABiunitIntervalGaussWeights30);
        break;
    case 31:
        setGaussPoints(IGABiunitIntervalGaussPoints31);
        setGaussWeights(IGABiunitIntervalGaussWeights31);
        break;
    case 32:
        setGaussPoints(IGABiunitIntervalGaussPoints32);
        setGaussWeights(IGABiunitIntervalGaussWeights32);
        break;
    case 33:
        setGaussPoints(IGABiunitIntervalGaussPoints33);
        setGaussWeights(IGABiunitIntervalGaussWeights33);
        break;
    case 34:
        setGaussPoints(IGABiunitIntervalGaussPoints34);
        setGaussWeights(IGABiunitIntervalGaussWeights34);
        break;
    case 35:
        setGaussPoints(IGABiunitIntervalGaussPoints35);
        setGaussWeights(IGABiunitIntervalGaussWeights35);
        break;
    case 36:
        setGaussPoints(IGABiunitIntervalGaussPoints36);
        setGaussWeights(IGABiunitIntervalGaussWeights36);
        break;
    case 37:
        setGaussPoints(IGABiunitIntervalGaussPoints37);
        setGaussWeights(IGABiunitIntervalGaussWeights37);
        break;
    case 38:
        setGaussPoints(IGABiunitIntervalGaussPoints38);
        setGaussWeights(IGABiunitIntervalGaussWeights38);
        break;
    case 39:
        setGaussPoints(IGABiunitIntervalGaussPoints39);
        setGaussWeights(IGABiunitIntervalGaussWeights39);
        break;
    case 40:
        setGaussPoints(IGABiunitIntervalGaussPoints40);
        setGaussWeights(IGABiunitIntervalGaussWeights40);
        break;
    case 41:
        setGaussPoints(IGABiunitIntervalGaussPoints41);
        setGaussWeights(IGABiunitIntervalGaussWeights41);
        break;
    case 42:
        setGaussPoints(IGABiunitIntervalGaussPoints42);
        setGaussWeights(IGABiunitIntervalGaussWeights42);
        break;
    case 43:
        setGaussPoints(IGABiunitIntervalGaussPoints43);
        setGaussWeights(IGABiunitIntervalGaussWeights43);
        break;
    case 44:
        setGaussPoints(IGABiunitIntervalGaussPoints44);
        setGaussWeights(IGABiunitIntervalGaussWeights44);
        break;
    case 45:
        setGaussPoints(IGABiunitIntervalGaussPoints45);
        setGaussWeights(IGABiunitIntervalGaussWeights45);
        break;
    case 46:
        setGaussPoints(IGABiunitIntervalGaussPoints46);
        setGaussWeights(IGABiunitIntervalGaussWeights46);
        break;
    case 47:
        setGaussPoints(IGABiunitIntervalGaussPoints47);
        setGaussWeights(IGABiunitIntervalGaussWeights47);
        break;
    case 48:
        setGaussPoints(IGABiunitIntervalGaussPoints48);
        setGaussWeights(IGABiunitIntervalGaussWeights48);
        break;
    case 49:
        setGaussPoints(IGABiunitIntervalGaussPoints49);
        setGaussWeights(IGABiunitIntervalGaussWeights49);
        break;
    case 50:
        setGaussPoints(IGABiunitIntervalGaussPoints50);
        setGaussWeights(IGABiunitIntervalGaussWeights50);
        break;
    default:
        assert(false);
    }
}

GaussQuadratureBiunitQuadrilateral::GaussQuadratureBiunitQuadrilateral(int _numGaussPoints) :
        GaussQuadrature(_numGaussPoints, 2) {
    switch (_numGaussPoints) {
    case 1:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints1);
        setGaussWeights(IGABiunitQuadrilateralWeights1);
        break;
    case 4:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints4);
        setGaussWeights(IGABiunitQuadrilateralWeights4);
        break;
    case 9:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints9);
        setGaussWeights(IGABiunitQuadrilateralWeights9);
        break;
    case 16:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints16);
        setGaussWeights(IGABiunitQuadrilateralWeights16);
        break;
    case 25:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints25);
        setGaussWeights(IGABiunitQuadrilateralWeights25);
        break;
    case 36:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints36);
        setGaussWeights(IGABiunitQuadrilateralWeights36);
        break;
    case 49:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints49);
        setGaussWeights(IGABiunitQuadrilateralWeights49);
        break;
    case 64:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints64);
        setGaussWeights(IGABiunitQuadrilateralWeights64);
        break;
    case 81:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints81);
        setGaussWeights(IGABiunitQuadrilateralWeights81);
        break;
    case 100:
        setGaussPoints(IGABiunitQuadrilateralGaussPoints100);
        setGaussWeights(IGABiunitQuadrilateralWeights100);
        break;
    default:
        assert(false);
    }
}

void computeLowOrderShapeFunctions(int _nNodes, const double *_coords, double *_shapeFuncs) {
    /*
     *  Returns the low order shape functions and their parametric derivatives as follows;
     * _shapeFuncs = [_shapeFuncs d(_shapeFuncs)/d(xi) d(_shapeFuncs)/d(eta)]
     */
    assert(_coords!=NULL);
    assert(_shapeFuncs!=NULL);
    if (_nNodes == 3) {
        // Shape function values
        _shapeFuncs[0] = 1 - _coords[0] - _coords[1];
        _shapeFuncs[1] = _coords[0];
        _shapeFuncs[2] = _coords[1];

        // Derivatives of the shape functions with respect to zeta1
        _shapeFuncs[3] = -1;
        _shapeFuncs[4] = 1;
        _shapeFuncs[5] = 0;

        // Derivatives of the shape functions with respect to zeta2
        _shapeFuncs[6] = -1;
        _shapeFuncs[7] = 0;
        _shapeFuncs[8] = 1;
    } else if (_nNodes == 4) {
        // Shape function values
        _shapeFuncs[0] = (1 - _coords[0]) / 2 * (1 - _coords[1]) / 2;
        _shapeFuncs[1] = (1 + _coords[0]) / 2 * (1 - _coords[1]) / 2;
        _shapeFuncs[2] = (1 + _coords[0]) / 2 * (1 + _coords[1]) / 2;
        _shapeFuncs[3] = (1 - _coords[0]) / 2 * (1 + _coords[1]) / 2;

        // Derivatives of the shape functions with respect to xi
        _shapeFuncs[4] = - (1 - _coords[1]) / 4;
        _shapeFuncs[5] = (1 - _coords[1]) / 4;
        _shapeFuncs[6] = (1 + _coords[1]) / 4;
        _shapeFuncs[7] = - (1 + _coords[1]) / 4;

        // Derivatives of the shape functions with respect to eta
        _shapeFuncs[8] = - (1 - _coords[0]) / 4;
        _shapeFuncs[9] = - (1 + _coords[0]) / 4;
        _shapeFuncs[10] = (1 + _coords[0]) / 4;
        _shapeFuncs[11] = (1 - _coords[0]) / 4;
    } else
        assert(false);
}

void computeLineareKombination(int _nNodes, int _nValue, const double *_values,
        const double *_shapeFuncs, double *_returnValue) {

    for (int i = 0; i < _nValue; i++) {
        _returnValue[i] = 0;
        for (int j = 0; j < _nNodes; j++)
            _returnValue[i] += _values[j * _nValue + i] * _shapeFuncs[j];
    }

}

/***********************************************************************************************
* \brief Compute the cross product between two vectors in the 3-D space
* \param[in/out] _product The product of vector1 and vector 2
* \param[in] _v1 The 1st vector
* \param[in] _v2 The 2nd vector
* \author Andreas Apostolatos
***********/
void crossProduct3D(double* _product, double* _u, double* _v) {
    // Check input
    assert(_product != NULL);
    assert(_u != NULL);
    assert(_v != NULL);

    // Compute the cross product using the permutation tensor
    _product[0] = _u[1] * _v[2] - _u[2] * _v[1];
    _product[1] = _u[2] * _v[0] - _u[0] * _v[2];
    _product[2] = _u[0] * _v[1] - _u[1] * _v[0];

    // No return value
    return;
}

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

double getGavitationalAcceleration(string text) {
    assert(text != "");
    stringstream ss(text);
    double omega;
    ss >> omega;
    return omega;
}

void computeDirectElemTable(vector<int>** directElemTable, int numElems, int numNodes, int* nodeIDs, int* numNodesPerElem, int* elemTable) {
    /*
     * Returns the direct element table
     */
    map<int, int> *nodesMap = new map<int, int>();

    for (int i = 0; i < numNodes; i++)
        nodesMap->insert(nodesMap->end(), pair<int, int>(nodeIDs[i], i));

    int count = 0;
    for (int i = 0; i < numElems; i++) {
        const int numNodesElem = numNodesPerElem[i];
        for (int j = 0; j < numNodesElem; j++) {
            directElemTable[i]->push_back(nodesMap->at(elemTable[count + j]));
        }
        count += numNodesElem;
    }
    delete nodesMap;
}

/***********************************************************************************************
 * \brief Returns the centrifufal forces for the two blade system
 * \param[in/out] forceVctCentrifugal The force vector related to the centrifugal forces
 * \param[in] omega The rotational acceleration
 * \param[in] density The density of the blades
 * \param[in] directElemTable Table relating the elements and the nodes using global numbering
 * \param[in] numElems Number of elements in the mesh
 * \param[in] numNodesPerElem Vector containing the number of nodes for each element
 * \param[in] nodeCoors Vector of doubles containing the coordinates of the nodes in the mesh
 * \author Andreas Apostolatos
 ***********/
void computeCentrifugalForces(double* forceVctCentrifugal, double density, double omega, double thickness, vector<int>** directElemTable, int numElems, int* numNodesPerElem, double* nodeCoors) {
    // Initialize variables
    const int numGPsOnTri = 3;
    const int numGPsOnQuad = 4;
    int noCoord = 3;
    double gaussWeight;
    double detJ;
    double cartesianCoord[noCoord];
    double dirVct[noCoord];
    double normalVct[noCoord];
    double baseVct1[noCoord];
    double baseVct2[noCoord];
    double center[3] = {0.0, 0.0, 0.0};

    // Loop over all elements
    for (int i = 0; i < numElems; i++) {
        // Get the number of the element nodes
        int numNodesElem = numNodesPerElem[i];
        if (numNodesElem != 3 && numNodesElem != 4)
            assert(false);
        double *elem = new double[numNodesElem * 3];

        // Initialize variables
        double basisFunctions[numNodesElem];
        double basisFunctionsAndDerivatives[numNodesElem * 3];

        // Get the coordinates of the nodes
        for (int j = 0; j < numNodesElem; j++) {
            int nodePos = directElemTable[i]->at(j); // position of the node
            for (int k = 0; k < noCoord; k++) {
                elem[j * 3 + k] = nodeCoors[nodePos * 3 + k];
            }
        }

        // Get the numbering of the nodes
        int posNodes[numNodesElem]; // the position of the nodes in the master element
        for (int j = 0; j < numNodesElem; j++)
            posNodes[j] = directElemTable[i]->at(j);

        // Get the Gauss quadrature
        GaussQuadrature* theGaussQuadrature;
        if (numNodesElem == 3) {
            int numGPs = numGPsOnTri;
            theGaussQuadrature = new GaussQuadratureTriangle(numGPs);
        } else {
            int numGPs = numGPsOnQuad;
            theGaussQuadrature = new GaussQuadratureBiunitQuadrilateral(numGPs);
        }

        // Loop over all Gauss points
        for (int j = 0; j < theGaussQuadrature->getNumGaussPoints(); j++) {
            // Get the Gauss point coordinates in the integration space
            const double *gaussPoint = theGaussQuadrature->getGaussPoint(j);

            // Get the Gauss point weight
            gaussWeight = theGaussQuadrature->getGaussWeight(j);

            // Compute the basis functions describing parametrically the integration domain
            computeLowOrderShapeFunctions(numNodesElem, gaussPoint, basisFunctionsAndDerivatives);
            for (int k = 0; k < numNodesElem; k++)
                basisFunctions[k] = basisFunctionsAndDerivatives[0*numNodesElem + k];

            // Compute the Jacobian determinant
            for (int l = 0; l < noCoord; l++) {
                baseVct1[l] = 0.0;
                baseVct2[l] = 0.0;
            }
            for (int k = 0; k < numNodesElem; k++)
                for (int l = 0; l < noCoord; l++) {
                    baseVct1[l] += basisFunctionsAndDerivatives[1*numNodesElem + k]*elem[k*3 + l];
                    baseVct2[l] += basisFunctionsAndDerivatives[2*numNodesElem + k]*elem[k*3 + l];
                }
            crossProduct3D(normalVct,baseVct1,baseVct2);
            detJ = 0.0;
            for (int k = 0; k < noCoord; k++)
                detJ += pow(normalVct[k],2.0);
            detJ = sqrt(detJ); // The Gauss points for a triangle already contain the x2 factor which should have been included here! Therefore it can be nicely grouped together with the quad

            // Get the Gauss point weight
            gaussWeight = theGaussQuadrature->getGaussWeight(j);

            // Compute the image of the Gauss point in the Cartesian space
            computeLineareKombination(numNodesElem, noCoord, elem, basisFunctions, cartesianCoord);

            // Compute the unit vector along which the centrifugal foce is applied
            for (int k = 0; k < noCoord; k++)
                dirVct[k] = cartesianCoord[k] - center[k];

            // Project the force to the x-axis
            dirVct[1] = 0.0;
            dirVct[2] = 0.0;

            // Loop over all nodes of the element and all Cartesian directions and add the contribution to the global centrifugal force vector
            for (int k = 0; k < numNodesElem; k++) {
                for (int l = 0; l < noCoord; l++) {
                    forceVctCentrifugal[posNodes[k] * 3 + l] += density*pow(omega,2.0)*dirVct[l]*basisFunctions[k]*thickness*detJ*gaussWeight;
                }
            }
        }
    }
}

/***********************************************************************************************
 * \brief Returns the gravitational forces for the two blade system
 * \param[in/out] forceVctCentrifugal The force vector related to the centrifugal forces
 * \param[in] density The density of the blades
 * \param[in] g The gravitational acceleration
 * \param[in] directElemTable Table relating the elements and the nodes using global numbering
 * \param[in] numElems Number of elements in the mesh
 * \param[in] numNodesPerElem Vector containing the number of nodes for each element
 * \param[in] nodeCoors Vector of doubles containing the coordinates of the nodes in the mesh
 * \author Andreas Apostolatos
 ***********/
void computeGravitationalForces(double* forceVctGavitational, double density, double g, double thickness, vector<int>** directElemTable, int numElems, int* numNodesPerElem, double* nodeCoors) {
    // Initialize variables
    const int numGPsOnTri = 3;
    const int numGPsOnQuad = 4; // 4
    int noCoord = 3;
    double gaussWeight;
    double detJ;
    double dirVct[3] = {0.0, 0.0, 1.0};
    double normalVct[noCoord];
    double baseVct1[noCoord];
    double baseVct2[noCoord];

    // Loop over all elements
    for (int i = 0; i < numElems; i++) {
        // Get the number of the element nodes
        int numNodesElem = numNodesPerElem[i];
        if (numNodesElem != 3 && numNodesElem != 4)
            assert(false);
        double *elem = new double[numNodesElem * 3];

        // Initialize variables
        double basisFunctions[numNodesElem * 3];

        // Get the coordinates of the nodes
        for (int j = 0; j < numNodesElem; j++) {
            int nodePos = directElemTable[i]->at(j); // position of the node
            for (int k = 0; k < noCoord; k++) {
                elem[j * 3 + k] = nodeCoors[nodePos * 3 + k];
            }
        }

        // Get the numbering of the nodes
        int posNodes[numNodesElem]; // the position of the nodes in the master element
        for (int j = 0; j < numNodesElem; j++)
            posNodes[j] = directElemTable[i]->at(j);

        // Get the Gauss quadrature
        GaussQuadrature* theGaussQuadrature;
        if (numNodesElem == 3) {
            int numGPs = numGPsOnTri;
            theGaussQuadrature = new GaussQuadratureTriangle(numGPs);
        } else {
            int numGPs = numGPsOnQuad;
            theGaussQuadrature = new GaussQuadratureBiunitQuadrilateral(numGPs);
        }

        // Loop over all Gauss points
        for (int j = 0; j < theGaussQuadrature->getNumGaussPoints(); j++) {
            // Get the Gauss point coordinates in the integration space
            const double *gaussPoint = theGaussQuadrature->getGaussPoint(j);

            // Get the Gauss point weight
            gaussWeight = theGaussQuadrature->getGaussWeight(j);

            // Compute the basis functions describing parametrically the integration domain
            computeLowOrderShapeFunctions(numNodesElem, gaussPoint, basisFunctions);

            // Compute the Jacobian determinant
            for (int l = 0; l < noCoord; l++) {
                baseVct1[l] = 0.0;
                baseVct2[l] = 0.0;
            }
            for (int k = 0; k < numNodesElem; k++)
                for (int l = 0; l < noCoord; l++) {
                    baseVct1[l] += basisFunctions[1*numNodesElem + k]*elem[k*3 + l];
                    baseVct2[l] += basisFunctions[2*numNodesElem + k]*elem[k*3 + l];
                }
            crossProduct3D(normalVct,baseVct1,baseVct2);
            detJ = 0.0;
            for (int k = 0; k < noCoord; k++)
                detJ += pow(normalVct[k],2.0);
            detJ = sqrt(detJ); // The Gauss points for a triangle already contain the x2 factor which should have been included here! Therefore it can be nicely grouped together with the quad

            // Get the Gauss point weight
            gaussWeight = theGaussQuadrature->getGaussWeight(j);

            // Loop over all nodes of the element and all Cartesian directions and add the contribution to the global centrifugal force vector
            for (int k = 0; k < numNodesElem; k++) {
                for (int l = 0; l < noCoord; l++) {
                    forceVctGavitational[posNodes[k] * 3 + l] += density*g*dirVct[l]*basisFunctions[0*numNodesElem + k]*thickness*detJ*gaussWeight;
                }
            }
        }
    }
}

void rotateForces(double phi, int numNodes, double *forces) {
    int noCoord = 3;
    phi = -phi;
    for (int i = 0; i < numNodes; i++) {
        double f_x = cos(phi) * forces[i * noCoord + 0] + sin(phi) * forces[i * noCoord + 2];
        double f_z = - sin(phi) * forces[i * noCoord + 0] + cos(phi) * forces[i * noCoord + 2];
        forces[i * noCoord + 0] = f_x;
        forces[i * noCoord + 2] = f_z;
    }
}

/*
 * NRELBlades reads GiD .msh and .res file, send result at each timestep to the server
 */
int main(int argc, char **argv) {
//    EMPIRE_API_Connect("NRELBlades.xml");
    if (argc != 2) {
        std::cout << "Provide a valid input file" << std::endl;
        exit(-1);
    } else
        EMPIRE_API_Connect(argv[1]);

    int numTimeSteps = getNumTimeSteps(EMPIRE_API_getUserDefinedText("numTimeSteps"));
    double density = getOmega(EMPIRE_API_getUserDefinedText("density"));
    double g = getGavitationalAcceleration(EMPIRE_API_getUserDefinedText("gravitationalAcceleration"));
    double omega = getOmega(EMPIRE_API_getUserDefinedText("omega"));
    double thickness = getOmega(EMPIRE_API_getUserDefinedText("thickness"));
    double timeInterval = getTimeInterval(EMPIRE_API_getUserDefinedText("timeInterval"));
    double noStartTimeStep = getTimeInterval(EMPIRE_API_getUserDefinedText("noStartTimeStep"));
    bool todoIterativeCoupling = doIterativeCoupling(EMPIRE_API_getUserDefinedText("couplingType"));

    // Read the mesh using a GiD mesh file
    int numNodes = -1;
    int numElems = -1;
    double *nodeCoors;
    int *nodeIDs;
    int *numNodesPerElem;
    int *elemTable;
    int *elemIDs;
    string meshFile = EMPIRE_API_getUserDefinedText("GiDMeshFile");
    GiDFileIO::readDotMsh(meshFile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
            elemTable, elemIDs);
    EMPIRE_API_sendMesh("", numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable);

    // Initialize auxiliary variables
    int noCoord = 3;
    double forceVctCentrifugal[noCoord * numNodes];
    double forceVctGavitational[noCoord * numNodes];
    double forceVct[noCoord * numNodes];

    // Compute the direct element table
    std::vector<int> ** directElemTable = new vector<int>*[numElems];
    for (int i = 0; i < numElems; i++)
        directElemTable[i] = new vector<int>;
    computeDirectElemTable(directElemTable, numElems, numNodes, nodeIDs, numNodesPerElem, elemTable);

    // Compute the centrifugal force
    for (int i = 0; i < 3 * numNodes; i++)
        forceVctCentrifugal[i] = 0.0;
    computeCentrifugalForces(forceVctCentrifugal, density, omega, thickness, directElemTable, numElems, numNodesPerElem, nodeCoors);

    // Compute the gravitational force
    for (int i = 0; i < 3 * numNodes; i++)
        forceVctGavitational[i] = 0.0;
    computeGravitationalForces(forceVctGavitational, density, g, thickness, directElemTable, numElems, numNodesPerElem, nodeCoors);

    // Delete pointers
    delete[] nodeCoors;
    delete[] nodeIDs;
    delete[] numNodesPerElem;
    delete[] elemTable;
    delete[] elemIDs;

    // Loop over all time steps
    for (int i = 0; i < numTimeSteps; i++) {
        double phi = omega * timeInterval * (double) (i + noStartTimeStep);
        std::cout << "Do time step " << i + 1 << " ..." << std::endl;
        do {
            // Print the angle of rotation
            std::cout << "Angle of rotation : " << phi << " ..." << std::endl;

            // Receive field from Empire
            EMPIRE_API_recvDataField("defaultField", numNodes * 3, forceVct);

            // Add the gravitational forces
            for (int j = 0; j < numNodes * 3; j++)
                forceVct[j] += forceVctGavitational[j];

            // Rotate the forces
            rotateForces(phi, numNodes, forceVct);

            // Add the centrifugal forces
            for (int j = 0; j < numNodes * 3; j++)
                forceVct[j] += forceVctCentrifugal[j];

            // Send field to Empire
            EMPIRE_API_sendDataField("defaultField", numNodes * 3, forceVct);

            // do not receive any data
            if (!todoIterativeCoupling) {
                break;
            }

        } while (!EMPIRE_API_recvConvergenceSignal());
    }

    EMPIRE_API_Disconnect();
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
