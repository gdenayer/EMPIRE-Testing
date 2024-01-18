#include "QuadMortarMapper.h"
#include <assert.h>
#include <iostream>
#include "MortarMath.h"
#ifdef INTEL_MKL
#include <mkl_pardiso.h>
#endif
#include <stdlib.h>

using namespace std;

/********//**
 * \brief Shape functions product on the same quadrilateral. (N1*N1, N1*N2, N1*N3, N1*N4, ...)
 ***********/
class QuadMortarMapper::ShapeFuncProductOnSameQuad: public Integrand {
private:
    /// 1st shape function index (0 or 1 or 2 or 3)
    int shapeFunc1;
    /// 2nd shape function index (0 or 1 or 2 or 3)
    int shapeFunc2;
    /// all shape function values (N1, N2, N3, N4) on all Gauss points (singleton)
    static double **GPShapeFunc;
    /// determinant of jaccobian on all gauss points
    double *detJs;

    /***********************************************************************************************
     * \brief compute GPShapeFunc.
     * \param[in] NUM number of Gauss points
     ***********/
    static void initGPShapeFunc(int NUM) {
        if (GPShapeFunc == NULL) {
            GPShapeFunc = new double*[NUM];
            for (int i = 0; i < NUM; i++)
                GPShapeFunc[i] = new double[4];
            for (int i = 0; i < NUM; i++) {
                const double *gaussPoint;
                gaussPoint = &(GaussQuadrature::getQuadGaussPoints(NUM)[i * 2]);
                for (int j = 0; j < 4; j++)
                    GPShapeFunc[i][j] = MortarMath::getQuadShapeFunc(j, gaussPoint);
            }
        }
    }

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] quad the quadrilateral where Gauss quadrature is computed
     ***********/
    ShapeFuncProductOnSameQuad(const double *quad) {
        const int NUM = 4; // number of gauss points
        setNumGaussPoints(NUM);
        initGPShapeFunc(NUM);
        detJs = new double[NUM];
        for (int i = 0; i < NUM; i++) {
            const double *gaussPoint;
            gaussPoint = &(GaussQuadrature::getQuadGaussPoints(NUM)[i * 2]);
            detJs[i] = MortarMath::computeQuadDetJ(quad, gaussPoint);
        }
    }

    /***********************************************************************************************
     * \brief Destructor
     ***********/
    virtual ~ShapeFuncProductOnSameQuad() {
        delete[] detJs;
    }

    /***********************************************************************************************
     * \brief Set to multiply which two shape functions (N? * N?).
     * \param[in] _shapeFunc1 index of the 1st shape function, could be 0, 1, 2 or 3.
     * \param[in] _shapeFunc2 index of the 2nd shape function, could be 0, 1, 2 or 3. If _shapeFunc2 < 0, integrate shapeFunc1.
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
        int i = gaussPointID;
        if (shapeFunc2 >= 0)
            return detJs[i] * GPShapeFunc[i][shapeFunc1] * GPShapeFunc[i][shapeFunc2];
        else
            return detJs[i] * GPShapeFunc[i][shapeFunc1];
    }
};

double **QuadMortarMapper::ShapeFuncProductOnSameQuad::GPShapeFunc;

/********//**
 * \brief Shape functions product on 2 quadrilaterals. (N_A1*N_B1, N_A1*N_B2, N_A1*N_B3, N_A1*N_B4, ...)
 ***********/
class QuadMortarMapper::ShapeFuncProductOn2Quads: public Integrand {
private:
    /// the 1st quadrilateral
    const double *quad1;
    /// the 2nd quadrilateral
    const double *quad2;
    /// when computing local coordinates, project to x-y or y-z or z-x plane (only use 2 coordinates)
    int planeToProject;
    /// a triangle of the clip (the clipped polygon is divided into triangles)
    const double *clipTriangle;
    /// the area of the clipTriangle
    double clipArea;
    /// the local coordinates in quad1 of the Gauss points got from clipTriangle
    double *gaussPointLocalCoorsInQuad1;
    /// the local coordinates in quad2 of the Gauss points got from clipTriangle
    double *gaussPointLocalCoorsInQuad2;
    /// the shape function values in quad1 of the Gauss points got from clipTriangle
    double *GPShapeFuncInQuad1;
    /// the shape function values in quad2 of the Gauss points got from clipTriangle
    double *GPShapeFuncInQuad2;
    /// index of the 1st shape function, could be 0, 1, 2 or 3.
    int shapeFunc1;
    /// index of the 2nd shape function, could be 0, 1, 2 or 3.
    int shapeFunc2;
    /// if there exists one Gauss point lying outside quad1 or quad2, set the integrand to 0
    bool gaussPointInsideQuad;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _quad1 the 1st quadrilateral
     * \param[in] _quad2 the 2nd quadrilateral
     * \param[in] _planeToProject when computing local coordinates, project to x-y or y-z or z-x plane (only use 2 coordinates)
     * \param[in] _clipTriangle a triangle of the clip (the clipped polygon is divided into triangles)
     ***********/
    ShapeFuncProductOn2Quads(const double *_quad1, const double *_quad2, int _planeToProject,
            const double *_clipTriangle) :
            quad1(_quad1), quad2(_quad2), planeToProject(_planeToProject), clipTriangle(
                    _clipTriangle) {
        // empty
        // Test Case
        const int NUM = 3; // the number of gauss points
        gaussPointLocalCoorsInQuad1 = new double[NUM * 2];
        gaussPointLocalCoorsInQuad2 = new double[NUM * 2];
        GPShapeFuncInQuad1 = new double[NUM * 4];
        GPShapeFuncInQuad2 = new double[NUM * 4];
        setNumGaussPoints(NUM);
        for (int i = 0; i < NUM; i++) {
            const double *gaussPoint = &(GaussQuadrature::getTriangleGaussPoints(NUM)[i * 3]);
            double globalCoor[3];
            MortarMath::computeGlobalCoorInTriangle(clipTriangle, gaussPoint, globalCoor);
            if (!MortarMath::computeLocalCoorInQuad(quad1, planeToProject, globalCoor,
                    &gaussPointLocalCoorsInQuad1[i * 2])) {
                gaussPointInsideQuad = false;
                return;
            }

            if (!MortarMath::computeLocalCoorInQuad(quad2, planeToProject, globalCoor,
                    &gaussPointLocalCoorsInQuad2[i * 2])) {
                gaussPointInsideQuad = false;
                return;
            }
        }
        for (int i = 0; i < NUM; i++) {
            double *gaussPointInQuad1 = &gaussPointLocalCoorsInQuad1[i * 2];
            double *gaussPointInQuad2 = &gaussPointLocalCoorsInQuad2[i * 2];
            for (int j = 0; j < 4; j++) {
                GPShapeFuncInQuad1[i * 4 + j] = MortarMath::getQuadShapeFunc(j, gaussPointInQuad1);
                GPShapeFuncInQuad2[i * 4 + j] = MortarMath::getQuadShapeFunc(j, gaussPointInQuad2);
            }
        }
        clipArea = MortarMath::computeArea(clipTriangle);
        gaussPointInsideQuad = true;
    }
    /***********************************************************************************************
     * \brief Destructor
     ***********/
    virtual ~ShapeFuncProductOn2Quads() {
        delete[] gaussPointLocalCoorsInQuad1;
        delete[] gaussPointLocalCoorsInQuad2;
        delete[] GPShapeFuncInQuad1;
        delete[] GPShapeFuncInQuad2;
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
        int i = gaussPointID;
        if (gaussPointInsideQuad)
            return clipArea * GPShapeFuncInQuad1[i * 4 + shapeFunc1]
                    * GPShapeFuncInQuad2[i * 4 + shapeFunc2];
        else
            return 0.0;
    }
};

QuadMortarMapper::QuadMortarMapper(int _slaveNumNodes, int _slaveNumElems, int _slaveNodesPerElem,
        double *_slaveNodeCoors, int *_slaveNodeNumbers, int *_slaveElemTable, int _masterNumNodes,
        int _masterNumElems, int _masterNodesPerElem, double *_masterNodeCoors,
        int *_masterNodeNumbers, int *_masterElemTable, bool _oppositeSurfaceNormal, bool _dual) :
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
    setN_P_E(4);
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

QuadMortarMapper::~QuadMortarMapper() {
    // empty
}

void QuadMortarMapper::computeSlaveElemNormals() {
    slaveElemNormals = new double[slaveNumElems * 3];
    for (int i = 0; i < slaveNumElems; i++) {
        double slaveQuad[12];
        int *pos = &slaveDirectElemTable[i * 4];
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 3; k++)
                slaveQuad[j * 3 + k] = slaveNodeCoors[pos[j] * 3 + k];
        MortarMath::computeQuadNormal(slaveQuad, true, &slaveElemNormals[i * 3]);
    }
}

Integrand *QuadMortarMapper::getShapeFuncProductOnSingleElem(const double *elem) {
    return new ShapeFuncProductOnSameQuad(elem);
}

Integrand *QuadMortarMapper::getShapeFuncProductOn2Elems(const double *elem1, const double *elem2,
        int planeToProject, const double *clipTriangle) {
    return new ShapeFuncProductOn2Quads(elem1, elem2, planeToProject, clipTriangle);
}
