#include "TriangleMortarMapper.h"
#include <iostream>
#include <stdlib.h>
#include "MortarMath.h"
#include <assert.h>
#include <iomanip>
#include <math.h>
#ifdef INTEL_MKL
#include <mkl_pardiso.h>
#endif
using namespace std;

/********//**
 * \brief Shape functions product on the same triangle. (N1*N1, N1*N2, N1*N3, N2*N2, N2*N3, N3*N3)
 ***********/
class TriangleMortarMapper::ShapeFuncProductOnSameTriangle: public Integrand {
private:
    /// area of the triangle
    double area;
    /// 1st shape function index (0 or 1 or 2)
    int shapeFunc1;
    /// 2nd shape function index (0 or 1 or 2)
    int shapeFunc2;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] triangle the triangle where Gauss quadrature is computed
     ***********/
    ShapeFuncProductOnSameTriangle(const double *triangle) {
        const int NUM = 3; // number of Gauss points
        area = MortarMath::computeArea(triangle);
        setNumGaussPoints(NUM);
    }
    /***********************************************************************************************
     * \brief Destructor
     ***********/
    virtual ~ShapeFuncProductOnSameTriangle() {
        // empty
    }
    /***********************************************************************************************
     * \brief Set to multiply which two shape functions (N? * N?).
     * \param[in] _shapeFunc1 index of the 1st shape function, could be 0, 1, 2.
     * \param[in] _shapeFunc2 index of the 2nd shape function, could be 0, 1, 2. If _shapeFunc2 < 0, integrate shapeFunc1.
     ***********/
    void setShapeFunctions(int _shapeFunc1, int _shapeFunc2) {
        shapeFunc1 = _shapeFunc1;
        shapeFunc2 = _shapeFunc2;
    }
    /***********************************************************************************************
     * \brief Evaluate the value of shape function product on a Gauss point
     * \param[in] gaussPointID id of the Gauss point.
     ***********/
    double evaluate(int gaussPointID) {
        // on an arbitrary point of a triangle, the i-th shape function value
        // is the same as the i-th isoparametric coordinate
        const int NUM = getNumGaussPoints();
        const double *gaussPoint = &(GaussQuadrature::getTriangleGaussPoints(NUM)[gaussPointID * 3]);
        if (shapeFunc2 >= 0)
            return area * gaussPoint[shapeFunc1] * gaussPoint[shapeFunc2];
        else
            return area * gaussPoint[shapeFunc1];
    }
};

/********//**
 * \brief Shape functions product on 2 triangles. (N_A1*N_B1, N_A1*N_B2, N_A1*N_B3, ...)
 ***********/
class TriangleMortarMapper::ShapeFuncProductOn2Triangles: public Integrand {
private:
    /// the 1st triangle
    const double *triangle1;
    /// the 2nd triangle
    const double *triangle2;
    /// when computing local coordinates, project to x-y or y-z or z-x plane (only use 2 coordinates)
    int planeToProject;
    /// a triangle of the clip (the clipped polygon is divided into triangles)
    const double *clipTriangle;
    /// the area of the clipTriangle
    double clipArea;
    /// the local coordinates in triangle1 of the Gauss points got from clipTriangle
    double *gaussPointLocalCoorsInTri1;
    /// the local coordinates in triangle2 of the Gauss points got from clipTriangle
    double *gaussPointLocalCoorsInTri2;
    /// index of the 1st shape function, could be 0, 1, 2.
    int shapeFunc1;
    /// index of the 2nd shape function, could be 0, 1, 2.
    int shapeFunc2;
    /// if there exists one Gauss point lying outside triangle1 or triangle2, set the integrand to 0
    bool gaussPointInsideTri;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _triangle1 the 1st triangle
     * \param[in] _triangle2 the 2nd triangle
     * \param[in] _planeToProject when computing local coordinates, project to x-y or y-z or z-x plane (only use 2 coordinates)
     * \param[in] _clipTriangle a triangle of the clip (the clipped polygon is divided into triangles)
     ***********/
    ShapeFuncProductOn2Triangles(const double *_triangle1, const double *_triangle2,
            int _planeToProject, const double *_clipTriangle) :
            triangle1(_triangle1), triangle2(_triangle2), planeToProject(_planeToProject), clipTriangle(
                    _clipTriangle) {
        // empty
        const int NUM = 3; // the number of gauss points
        gaussPointLocalCoorsInTri1 = new double[NUM * 3];
        gaussPointLocalCoorsInTri2 = new double[NUM * 3];
        setNumGaussPoints(NUM);
        for (int i = 0; i < NUM; i++) {
            const double *gaussPoint = &(GaussQuadrature::getTriangleGaussPoints(NUM)[i * 3]);
            double globalCoor[3];
            MortarMath::computeGlobalCoorInTriangle(clipTriangle, gaussPoint, globalCoor);
            if (!MortarMath::computeLocalCoorInTriangle(triangle1, planeToProject, globalCoor,
                    &gaussPointLocalCoorsInTri1[i * 3])) {
                gaussPointInsideTri = false;
                return;
            }
            if (!MortarMath::computeLocalCoorInTriangle(triangle2, planeToProject, globalCoor,
                    &gaussPointLocalCoorsInTri2[i * 3])) {
                gaussPointInsideTri = false;
                return;
            }
        }
        clipArea = MortarMath::computeArea(clipTriangle);
        gaussPointInsideTri = true;
    }

    /***********************************************************************************************
     * \brief Destructor
     ***********/
    virtual ~ShapeFuncProductOn2Triangles() {
        delete[] gaussPointLocalCoorsInTri1;
        delete[] gaussPointLocalCoorsInTri2;
    }
    /***********************************************************************************************
     * \brief Set to multiply which two shape functions (N? * N?).
     * \param[in] _shapeFunc1 index of the 1st shape function, could be 0, 1, 2.
     * \param[in] _shapeFunc2 index of the 2nd shape function, could be 0, 1, 2.
     ***********/
    void setShapeFunctions(int _shapeFunc1, int _shapeFunc2) {
        shapeFunc1 = _shapeFunc1;
        shapeFunc2 = _shapeFunc2;
    }
    /***********************************************************************************************
     * \brief Evaluate the value of shape function product on a Gauss point
     * \param[in] gaussPointID id of the Gauss point.
     ***********/
    double evaluate(int gaussPointID) {
        // on an arbitrary point of a triangle, the i-th shape function value
        // is the same as the i-th isoparametric coordinate
        if (gaussPointInsideTri)
            return clipArea * gaussPointLocalCoorsInTri1[gaussPointID * 3 + shapeFunc1]
                    * gaussPointLocalCoorsInTri2[gaussPointID * 3 + shapeFunc2];
        else
            return 0.0;
    }
};

TriangleMortarMapper::TriangleMortarMapper(int _slaveNumNodes, int _slaveNumElems,
        int _slaveNodesPerElem, double *_slaveNodeCoors, int *_slaveNodeNumbers,
        int *_slaveElemTable, int _masterNumNodes, int _masterNumElems, int _masterNodesPerElem,
        double *_masterNodeCoors, int *_masterNodeNumbers, int *_masterElemTable,
        bool _oppositeSurfaceNormal, bool _dual, bool symmetric) :
        AbstractMortarMapper(_slaveNumNodes, _slaveNumElems, _slaveNodesPerElem, _slaveNodeCoors,
                _slaveNodeNumbers, _slaveElemTable, _masterNumNodes, _masterNumElems,
                _masterNodesPerElem, _masterNodeCoors, _masterNodeNumbers, _masterElemTable,
                _oppositeSurfaceNormal, _dual) {

    /* A - the slave, B - the master
     * C_BB and C_BA are the matrixes to be computed, and say the DOF vector
     * from A is W_A, the DOF vector from B is W_B, the mapping equation is:
     * C_BB * W_B = C_BA * W_A
     * the constructor will set up C_BB and C_BA
     */

    // 0. convert the quadrilateral mesh to triangular mesh
    convertToTriangularMesh(symmetric);
    setN_P_E(3);

    // 1. initialize assistant tables
    initTables();
    initANNTree();

    // 2. compute C_BB
    //cout<<"doing c_bb"<<endl;
    computeC_BB();

    // 3. compute C_BA
    //cout<<"doing c_ba"<<endl;
    computeC_BA();

    deleteANNTree();
    deleteTables();

    // 4. use pardiso to do factorization on C_BB
    if (!dual)
        initPardiso();

    // 5. check NULL pointers
    checkNullPointers();
}

TriangleMortarMapper::~TriangleMortarMapper() {
    // empty
}

void TriangleMortarMapper::convertToTriangularMesh(bool symmetric) {
    if (slaveNodesPerElem != 3 && masterNodesPerElem != 3 && slaveNodesPerElem != 4
            && masterNodesPerElem != 4) {
        cerr
                << "Error in TriangleMortarMapper: wrong mesh!!! Only allow triangles and quadrilaterals!!!"
                << endl;
        exit(EXIT_FAILURE);
    }

    if (symmetric) {
        factorC_BA = 1.0;
        factorC_BB = 1.0;
        if (masterNodesPerElem == 4) {
            int *elemTable = new int[masterNumElems * 3 * 4];
            for (int i = 0; i < masterNumElems; i++) {
                for (int j = 0; j < 3; j++) {
                    elemTable[4 * i * 3 + j] = masterElemTable[i * 4 + j];
                    elemTable[(4 * i + 1) * 3 + j] = masterElemTable[i * 4 + (j + 2) % 4];
                    elemTable[(4 * i + 2) * 3 + j] = masterElemTable[i * 4 + j + 1];
                    elemTable[(4 * i + 3) * 3 + j] = masterElemTable[i * 4 + (j + 3) % 4];
                }
            }
            masterElemTable = elemTable;
            masterNodesPerElem = 3;
            masterNumElems *= 4;
            factorC_BB = 2.0;
        }
        if (slaveNodesPerElem == 4) {
            int *elemTable = new int[slaveNumElems * 3 * 4];
            for (int i = 0; i < slaveNumElems; i++) {
                for (int j = 0; j < 3; j++) {
                    elemTable[4 * i * 3 + j] = slaveElemTable[i * 4 + j];
                    elemTable[(4 * i + 1) * 3 + j] = slaveElemTable[i * 4 + (j + 2) % 4];
                    elemTable[(4 * i + 2) * 3 + j] = slaveElemTable[i * 4 + j + 1];
                    elemTable[(4 * i + 3) * 3 + j] = slaveElemTable[i * 4 + (j + 3) % 4];
                }
            }
            slaveElemTable = elemTable;
            slaveNodesPerElem = 3;
            slaveNumElems *= 4;
            factorC_BA = 2.0;
        }
        factorC_BA *= factorC_BB;
    } else {
        factorC_BA = 1.0;
        factorC_BB = 1.0;
        if (slaveNodesPerElem == 4) {
            int *elemTable = new int[slaveNumElems * 3 * 2];
            for (int i = 0; i < slaveNumElems; i++) {
                for (int j = 0; j < 3; j++) {
                    elemTable[2 * i * 3 + j] = slaveElemTable[i * 4 + j];
                    elemTable[(2 * i + 1) * 3 + j] = slaveElemTable[i * 4 + (j + 2) % 4];
                }
            }
            slaveElemTable = elemTable;
            slaveNodesPerElem = 3;
            slaveNumElems *= 2;
        }
        if (masterNodesPerElem == 4) {
            int *elemTable = new int[masterNumElems * 3 * 2];
            for (int i = 0; i < masterNumElems; i++) {
                for (int j = 0; j < 3; j++) {
                    elemTable[2 * i * 3 + j] = masterElemTable[i * 4 + j];
                    elemTable[(2 * i + 1) * 3 + j] = masterElemTable[i * 4 + (j + 2) % 4];
                }
            }
            masterElemTable = elemTable;
            masterNodesPerElem = 3;
            masterNumElems *= 2;
        }
    }
}

void TriangleMortarMapper::computeSlaveElemNormals() {
    slaveElemNormals = new double[slaveNumElems * 3];
    for (int i = 0; i < slaveNumElems; i++) {
        double slaveTriangle[9];
        int *pos = &slaveDirectElemTable[i * 3];
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                slaveTriangle[j * 3 + k] = slaveNodeCoors[pos[j] * 3 + k];
        MortarMath::computeTriangleNormal(slaveTriangle, true, &slaveElemNormals[i * 3]);
    }
}

Integrand *TriangleMortarMapper::getShapeFuncProductOnSingleElem(const double *elem) {
    return new ShapeFuncProductOnSameTriangle(elem);
}

Integrand *TriangleMortarMapper::getShapeFuncProductOn2Elems(const double *elem1,
        const double *elem2, int planeToProject, const double *clipTriangle) {
    return new ShapeFuncProductOn2Triangles(elem1, elem2, planeToProject, clipTriangle);
}
