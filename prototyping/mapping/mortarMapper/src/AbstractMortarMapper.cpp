#include "AbstractMortarMapper.h"

#ifdef INTEL_MKL
#include <mkl_pardiso.h>
#include <mkl_spblas.h>
#include <mkl_lapacke.h>
#include <mkl_service.h>
#define LAPACK
#endif

#ifdef LAPACK_GCC
#include "lapacke.h"
#define LAPACK
#endif

#include <iostream>
#include <stdlib.h>
#include <set>
#include <assert.h>
#include <omp.h>

#include "MortarMath.h"
#include <math.h>

#ifdef FLANN
#include <flann/flann.hpp>
#endif
#ifdef ANN
#include "ANN/ANN.h"
#endif

using namespace std;

// set default value for MKL threads
int AbstractMortarMapper::mklSetNumThreads = 1;
// set default value for mapper threads
int AbstractMortarMapper::mapperSetNumThreads = 1;

AbstractMortarMapper::PolygonPointClipper::PolygonPointClipper(const double *_polygon, int _size,
        int _planeToProject) :
        polygon(_polygon), size(_size), planeToProject(_planeToProject) {

    // numbering the edge and point in the following way:
    // The edge i is the edge between point i and i+1
    //const double EPS = 1E-5; // decide whether to use y=f(x) or x=f(y)
    insideFlag = new bool[size];
    edgeSlopes = new double[size];
    reverseXY = new bool[size];

    for (int i = 0; i < size; i++) {
        // i is the ID of the edge
        double x1, y1, x2, y2, x3, y3;
        int x_pos, y_pos; // x, y index of all points in a triangle
        x_pos = (planeToProject + 1) % 3;
        y_pos = (planeToProject + 2) % 3;
        int p1_pos, p2_pos, p3_pos; // index of points of the edge in the triangle
        p1_pos = i;
        p2_pos = (i + 1) % size;
        p3_pos = (i + 2) % size;
        x1 = polygon[p1_pos * 3 + x_pos];
        y1 = polygon[p1_pos * 3 + y_pos];
        x2 = polygon[p2_pos * 3 + x_pos];
        y2 = polygon[p2_pos * 3 + y_pos];
        x3 = polygon[p3_pos * 3 + x_pos];
        y3 = polygon[p3_pos * 3 + y_pos];
        if (fabs(x2 - x1) < fabs(y2 - y1)) {
            reverseXY[i] = true;
            edgeSlopes[i] = (x2 - x1) / (y2 - y1);
            insideFlag[i] = (x3 - x1) > edgeSlopes[i] * (y3 - y1);
        } else {
            reverseXY[i] = false;
            edgeSlopes[i] = (y2 - y1) / (x2 - x1);
            insideFlag[i] = (y3 - y1) > edgeSlopes[i] * (x3 - x1);
        }
    }
}

AbstractMortarMapper::PolygonPointClipper::~PolygonPointClipper() {
    delete[] insideFlag;
    delete[] edgeSlopes;
    delete[] reverseXY;
}

bool AbstractMortarMapper::PolygonPointClipper::inside(int edgeID, double *point) {
    int x_pos, y_pos; // x, y pos of all points in a triangle
    x_pos = (planeToProject + 1) % 3;
    y_pos = (planeToProject + 2) % 3;
    int p1_pos = edgeID;
    double x1 = polygon[p1_pos * 3 + x_pos];
    double y1 = polygon[p1_pos * 3 + y_pos];
    double x4 = point[x_pos];
    double y4 = point[y_pos];
    if (reverseXY[edgeID])
        return ((x4 - x1) > edgeSlopes[edgeID] * (y4 - y1)) == insideFlag[edgeID];
    else
        return ((y4 - y1) > edgeSlopes[edgeID] * (x4 - x1)) == insideFlag[edgeID];
}

AbstractMortarMapper::AbstractMortarMapper(int _slaveNumNodes, int _slaveNumElems,
        int _slaveNodesPerElem, double *_slaveNodeCoors, int *_slaveNodeNumbers,
        int *_slaveElemTable, int _masterNumNodes, int _masterNumElems, int _masterNodesPerElem,
        double *_masterNodeCoors, int *_masterNodeNumbers, int *_masterElemTable,
        bool _oppositeSurfaceNormal, bool _dual) :
        slaveNumNodes(_slaveNumNodes), slaveNumElems(_slaveNumElems), slaveNodesPerElem(
                _slaveNodesPerElem), slaveNodeCoors(_slaveNodeCoors), slaveNodeNumbers(
                _slaveNodeNumbers), slaveElemTable(_slaveElemTable), masterNumNodes(
                _masterNumNodes), masterNumElems(_masterNumElems), masterNodesPerElem(
                _masterNodesPerElem), masterNodeCoors(_masterNodeCoors), masterNodeNumbers(
                _masterNodeNumbers), masterElemTable(_masterElemTable), oppositeSurfaceNormal(
                _oppositeSurfaceNormal), dual(_dual) {
    // a could-be NULL pointer must be initialized to NULL, otherwise there could be segmentation fault when delete it
    C_BB_A = NULL;
    C_BB_IA = NULL;
    C_BB_JA = NULL;
    C_BB_A_DUAL = NULL;
    C_BA_A = NULL;
    C_BA_IA = NULL;
    C_BA_JA = NULL;
    C_BA_A_DUAL = NULL;
    // set factorC_BB and factorC_BA
    factorC_BB = 1.0;
    factorC_BA = 1.0;
    // check whether the necessary libraries are there
#ifndef LAPACK
    cerr << endl;
    cerr
            << "AbstractMortarMapper::AbstractMortarMapper: No lapack library is found, mortar mapper cannot be used!"
            << endl;
    cerr << endl;
    exit(EXIT_FAILURE);
#endif
#ifndef INTEL_MKL
    if (!dual) {
        cerr << endl;
        cerr
                << "AbstractMortarMapper::AbstractMortarMapper: No pardiso library is found, standard mortar mapper cannot be used!"
                << endl << "\t Try dual mortar mapper!" << endl;
        exit(EXIT_FAILURE);
    }
#endif
}

AbstractMortarMapper::~AbstractMortarMapper() {
    /*delete[] slaveNodeCoors;
     delete[] slaveNodeNumbers;
     delete[] slaveElemTable;
     delete[] masterNodeCoors;
     delete[] masterNodeNumbers;
     delete[] masterElemTable;*/ // should not be deleted here
    if (!dual)
        deletePardiso();
    delete[] C_BB_A;
    delete[] C_BB_IA;
    delete[] C_BB_JA;
    delete[] C_BB_A_DUAL;
    delete[] C_BA_A;
    delete[] C_BA_IA;
    delete[] C_BA_JA;
    delete[] C_BA_A_DUAL;
}

void AbstractMortarMapper::consistentMapping(double *slaveDOF, double *masterDOF) {
    // 1. matrix vector product (W_tmp = C_BA * W_A)
    int m = masterNumNodes; // number of rows of C_BA
    int n = slaveNumNodes; // number of columns of C_BA
    char noTrans = 'N';
    char descra[] = "G00F"; // general matrix, indexing from 1
    double alpha = 1.0;
    double beta = 0.0;
#ifdef INTEL_MKL
    if (!dual) {
        mkl_set_num_threads(mklSetNumThreads);
        mkl_dcsrmv(&noTrans, &m, &n, &alpha, descra, C_BA_A, C_BA_JA, C_BA_IA, &C_BA_IA[1], slaveDOF,
                &beta, masterDOF);
    }
    else {
        mkl_set_num_threads(mklSetNumThreads);
        mkl_dcsrmv(&noTrans, &m, &n, &alpha, descra, C_BA_A_DUAL, C_BA_JA, C_BA_IA, &C_BA_IA[1], slaveDOF,
                &beta, masterDOF);
    }
#else
    if (!dual)
        MortarMath::dcsrmv(noTrans, m, n, C_BA_A, C_BA_JA, C_BA_IA, slaveDOF, masterDOF);
    else
        MortarMath::dcsrmv(noTrans, m, n, C_BA_A_DUAL, C_BA_JA, C_BA_IA, slaveDOF, masterDOF);
#endif

    // 2. solve C_BB * W_B = W_tmp
    if (!dual) {
        // pardiso forward and backward substitution
        int phase = 33; // forward and backward substitution
        int idum; // integer dummy
        double *ddum = new double[masterNumNodes]; // dummy but the memory is asked for
        int error = 0;
        iparm[5] = 1; // write solution to b
#ifdef INTEL_MKL
                mkl_set_num_threads(mklSetNumThreads);
                pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, C_BB_A, C_BB_IA, C_BB_JA, &idum, &nrhs,
                        iparm, &msglvl, masterDOF, ddum, &error);
#endif
        if (error != 0) {
            cerr << "Error in MortarMapper: pardiso solving failed!" << error << endl;
            exit(EXIT_FAILURE);
        }
        delete[] ddum;
    } else {
        for (int i = 0; i < masterNumNodes; i++)
            masterDOF[i] /= C_BB_A_DUAL[i];
    }
}

void AbstractMortarMapper::consistentMappingTraction2Force(double *slaveDOF, double *masterDOF) {
    // 1. matrix vector product (C_BB * t_B = C_BA * t_A ====> f_B = C_BA * t_A)
    int m = masterNumNodes; // number of rows of C_BA
    int n = slaveNumNodes; // number of columns of C_BA
    char noTrans = 'N';
    char descra[] = "G00F"; // general matrix, indexing from 1
    double alpha = 1.0;
    double beta = 0.0;
#ifdef INTEL_MKL
    mkl_set_num_threads(mklSetNumThreads);
    mkl_dcsrmv(&noTrans, &m, &n, &alpha, descra, C_BA_A, C_BA_JA, C_BA_IA, &C_BA_IA[1], slaveDOF,
            &beta, masterDOF);
#else
    MortarMath::dcsrmv(noTrans, m, n, C_BA_A, C_BA_JA, C_BA_IA, slaveDOF, masterDOF);
#endif
}

void AbstractMortarMapper::conservativeMapping(double *masterDOF, double *slaveDOF) {
    /*
     * consistent mapping:
     * C_BB * W_B = C_BA * W_A
     * since conservative mapping -- W_B^T * F_B = W_A^T * F_A
     * => F_A = C_BA^T * C_BB^(-1) * F_B
     * F_A --- slaveDOF
     * F_B --- masterDOF
     */
    // 1. solve C_BB * F_tmp = F_B
    if (!dual) {
        // pardiso forward and backward substitution
        int phase = 33; // forward and backward substitution
        int idum; // integer dummy
        double *ddum = new double[masterNumNodes]; // dummy but the memory is asked for
        int error = 0;
        iparm[5] = 1; // write solution to x
#ifdef INTEL_MKL
                mkl_set_num_threads(mklSetNumThreads);
                pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, C_BB_A, C_BB_IA, C_BB_JA, &idum, &nrhs,
                        iparm, &msglvl, masterDOF, ddum, &error);
#endif
        if (error != 0) {
            cerr << "Error in MortarMapper: pardiso solving failed!" << error << endl;
            exit(EXIT_FAILURE);
        }
        delete ddum;
    } else {
        for (int i = 0; i < masterNumNodes; i++)
            masterDOF[i] /= C_BB_A_DUAL[i];
    }

    // 2. matrix vector product (F_A = C_BA^T * F_tmp)
    int m = masterNumNodes; // number of rows of C_BA
    int n = slaveNumNodes; // number of columns of C_BA
    char trans = 'T';
    char descra[] = "G00F"; // general matrix, indexing from 1
    double alpha = 1.0;
    double beta = 0.0;
#ifdef INTEL_MKL
    if (!dual) {
        mkl_set_num_threads(mklSetNumThreads);
        mkl_dcsrmv(&trans, &m, &n, &alpha, descra, C_BA_A, C_BA_JA, C_BA_IA, &C_BA_IA[1], masterDOF,
                &beta, slaveDOF);
    }
    else {
        mkl_set_num_threads(mklSetNumThreads);
        mkl_dcsrmv(&trans, &m, &n, &alpha, descra, C_BA_A_DUAL, C_BA_JA, C_BA_IA, &C_BA_IA[1], masterDOF,
                &beta, slaveDOF);
    }
#else
    if (!dual)
        MortarMath::dcsrmv(trans, m, n, C_BA_A, C_BA_JA, C_BA_IA, masterDOF, slaveDOF);
    else
        MortarMath::dcsrmv(trans, m, n, C_BA_A_DUAL, C_BA_JA, C_BA_IA, masterDOF, slaveDOF);
#endif
}

void AbstractMortarMapper::conservativeMappingTraction2Force(double *masterDOF, double *slaveDOF) {
    /*
     * consistent mapping:
     * C_BB * W_B = C_BA * W_A
     * since conservative mapping -- W_B^T * F_B = W_A^T * F_A
     * => F_A = C_BA^T * C_BB^(-1) * F_B
     * F_A --- slaveDOF
     * F_B --- masterDOF
     * =========================================================
     * what is more, F_B = C_BB * t_B
     * => F_A = C_BA^T * t_B
     * REMARK: only C_BA is allowed, instead of C_BA_DUAL, this is why we store C_BA in case of dual
     */

    // 1. matrix vector product (F_A = C_BA^T * p_B)
    int m = masterNumNodes; // number of rows of C_BA
    int n = slaveNumNodes; // number of columns of C_BA
    char trans = 'T';
    char descra[] = "G00F"; // general matrix, indexing from 1
    double alpha = 1.0;
    double beta = 0.0;
#ifdef INTEL_MKL
    mkl_set_num_threads(mklSetNumThreads);
    mkl_dcsrmv(&trans, &m, &n, &alpha, descra, C_BA_A, C_BA_JA, C_BA_IA, &C_BA_IA[1], masterDOF,
            &beta, slaveDOF);
#else
    MortarMath::dcsrmv(trans, m, n, C_BA_A, C_BA_JA, C_BA_IA, masterDOF, slaveDOF);
#endif
}

void AbstractMortarMapper::computeC_BB() {
    // 1. compute the sparsity map
    // sparsity map has the information of a, ia, ja in a CSR formated matrix
    map<int, double> **sparsityMapC_BB = NULL; // a "could be" NULL pointer is better to be init to NULL
    if (!dual) {
        sparsityMapC_BB = new map<int, double>*[masterNumNodes];
        for (int i = 0; i < masterNumNodes; i++)
            sparsityMapC_BB[i] = new map<int, double>;
    } else {
        C_BB_A_DUAL = new double[masterNumNodes];
        for (int i = 0; i < masterNumNodes; i++)
            C_BB_A_DUAL[i] = 0.0;
    }

    for (int i = 0; i < masterNumElems; i++) {
        double elem[N_P_E * 3]; // this element
        int pos[N_P_E]; // the position of the node in the nodeCoors
        for (int j = 0; j < N_P_E; j++) {
            pos[j] = masterDirectElemTable[i * N_P_E + j];
            for (int k = 0; k < 3; k++)
                elem[j * 3 + k] = masterNodeCoors[pos[j] * 3 + k];
        }

        if (N_P_E == 4) { // replace the master element by the projection of it on its "element plane"
            double masterElemNormal[3];
            MortarMath::computeQuadNormal(elem, true, masterElemNormal);
            double masterQuadCenter[3];
            MortarMath::getQuadCenter(elem, masterQuadCenter);
            double masterQuadPrj[12];
            MortarMath::projectToElemPlane(masterQuadCenter, masterElemNormal, elem, 4,
                    masterQuadPrj);
            for (int i = 0; i < 12; i++)
                elem[i] = masterQuadPrj[i];
        }

        Integrand *integrand = getShapeFuncProductOnSingleElem(elem);
        // make use of the symmetry
        if (!dual) {
            for (int j = 0; j < N_P_E; j++) {
                for (int k = j; k < N_P_E; k++) {
                    integrand->setShapeFunctions(j, k);
                    double tmp = 0.0;
                    if (N_P_E == 3)
                        tmp = GaussQuadrature::gaussQuadratureOnTriangle(integrand);
                    else if (N_P_E == 4)
                        tmp = GaussQuadrature::gaussQuadratureOnQuad(integrand);

                    int smaller, larger;
                    if (pos[j] > pos[k]) {
                        larger = pos[j];
                        smaller = pos[k];
                    } else {
                        larger = pos[k];
                        smaller = pos[j];
                    }
                    bool inserted =
                            sparsityMapC_BB[smaller]->insert(pair<int, double>(larger, tmp)).second; // *.second is a bool indicating inserted or not
                    if (!inserted)
                        (sparsityMapC_BB[smaller]->at(larger)) += tmp; // at() returns a reference, so using "+=" is correct
                }
            }
        } else {
            for (int j = 0; j < N_P_E; j++) {
                integrand->setShapeFunctions(j, -1);
                double tmp = 0.0;
                if (N_P_E == 4)
                    tmp = GaussQuadrature::gaussQuadratureOnQuad(integrand);
                else if (N_P_E == 3)
                    tmp = GaussQuadrature::gaussQuadratureOnTriangle(integrand);
                C_BB_A_DUAL[pos[j]] += tmp;
            }
        }
        delete integrand;
    }

    // 2. according to sparsity map, compute a, ia, ja of CSR formated C_BB
    if (!dual) {
        C_BB_IA = new int[masterNumNodes + 1];
        C_BB_IA[0] = 1; // numbering from 1
        int nnz = 1; // number of non-zero entries
        for (int i = 0; i < masterNumNodes; i++) {
            nnz += sparsityMapC_BB[i]->size();
            C_BB_IA[i + 1] = nnz;
        }
        nnz--;
        C_BB_JA = new int[nnz];
        C_BB_A = new double[nnz];
        int count = 0;
        for (int i = 0; i < masterNumNodes; i++) {
            for (map<int, double>::iterator it = sparsityMapC_BB[i]->begin();
                    it != sparsityMapC_BB[i]->end(); it++) {
                C_BB_JA[count] = it->first + 1; // numbering from 1
                C_BB_A[count] = it->second / factorC_BB;
                count++;
            }
        }assert(nnz == count);

        // delete spasity map
        for (int i = 0; i < masterNumNodes; i++)
            delete sparsityMapC_BB[i];
        delete[] sparsityMapC_BB;

        /*cout << endl << "=====C_BB_IA======" << endl;
         for (int i = 0; i < masterNumNodes + 1; i++)
         cout << C_BB_IA[i] << "   ";
         cout << endl << "=====C_BB_JA======" << endl;
         for (int i = 0; i < nnz; i++)
         cout << C_BB_JA[i] << "   ";
         cout << endl << "=====C_BB_A======" << endl;
         for (int i = 0; i < nnz; i++)
         cout << C_BB_A[i] << "   ";

         double aaa[4][4];
         int pp = 0;
         for (int i=0; i<4; i++)
         for (int j=0; j<4; j++)
         if (i <= j) {
         aaa[i][j] = C_BB_A[pp];
         aaa[j][i] = C_BB_A[pp];
         pp++;
         }

         double ttt[4] = {0,0,0,0};
         for (int i=0; i<4; i++)
         for (int j=0; j<4; j++)
         ttt[i]+=aaa[i][j];
         cout << endl << "++++++++ " << endl;
         for (int j=0; j<4; j++)
         cout << ttt[j] << endl;*/
    } else {
        for (int i = 0; i < masterNumNodes; i++)
            C_BB_A_DUAL[i] /= factorC_BB;
    }
}

void AbstractMortarMapper::computeC_BA() {
    // 1. construct a nearest neighbor searching tree by the slave nodes
    // done outside

    // 2. initialize sparsity map, which has the information of a, ia, ja in a CSR formated matrix
    map<int, double> **sparsityMapC_BA = new map<int, double>*[masterNumNodes];
    for (int i = 0; i < masterNumNodes; i++)
        sparsityMapC_BA[i] = new map<int, double>;
    map<int, double> **sparsityMapC_BA_DUAL = NULL;
    if (dual) {
        sparsityMapC_BA_DUAL = new map<int, double>*[masterNumNodes];
        for (int i = 0; i < masterNumNodes; i++)
            sparsityMapC_BA_DUAL[i] = new map<int, double>;
    }

    // 3. compute entries in the sparsity map by looping over the master elements
#ifdef FLANN
#pragma omp parallel num_threads(mapperSetNumThreads)
#endif
    {
#ifdef FLANN
        //EMPIRE::AuxiliaryFunctions::report_num_threads(2);
#pragma omp for
#endif
        for (int i = 0; i < masterNumElems; i++) {
            // 3.1 compute the searching radius
            double masterElem[N_P_E * 3];
            getElemCoor(i, AbstractMortarMapper::MASTER, masterElem);
            double radiusSqr = computeSearchRadiusSquare(masterElem);
            // if (!(i%1000)) cout << endl << i << "RADIUS: " << radius << endl;

            // 3.2 find all candidates which may overlap the master element
            double masterElemNormal[3];
            if (N_P_E == 3) {
                MortarMath::computeTriangleNormal(masterElem, true, masterElemNormal);
            } else {
                MortarMath::computeQuadNormal(masterElem, true, masterElemNormal);
            }
            set<int> neighborElems;
            findCandidates(masterElem, masterElemNormal, radiusSqr, neighborElems);

            map<int, double*> projections; // the projections of all neighboring nodes
            if (N_P_E == 3)
                projectToElemPlane(masterElem, masterElemNormal, neighborElems, projections);
            if (N_P_E == 4) { // replace the master element by the projection of it on its "element plane"
                double masterQuadCenter[3];
                MortarMath::getQuadCenter(masterElem, masterQuadCenter);
                double masterQuadPrj[12];
                MortarMath::projectToElemPlane(masterQuadCenter, masterElemNormal, masterElem, 4,
                        masterQuadPrj);
                for (int j = 0; j < 12; j++)
                    masterElem[j] = masterQuadPrj[j];
                projectToElemPlane(masterQuadCenter, masterElemNormal, neighborElems, projections);
            }

            // if (!(i%1000)) cout << "candidates: " << neighborElems.size() << endl;

            int posMasterNodes[N_P_E]; // the position of the nodes in the master quad
            for (int j = 0; j < N_P_E; j++)
                posMasterNodes[j] = masterDirectElemTable[i * N_P_E + j];

            // 3.3 create the point clipper
            int planeToProject = MortarMath::projectToPlane(masterElemNormal);
            AbstractMortarMapper::PolygonPointClipper clipper(masterElem, N_P_E, planeToProject);

            // 3.3.5 if dual, compute the coefficient matrix here
            double coeffMatrix[16];
            if (dual)
                computeDualCoeffMatrix(masterElem, coeffMatrix);

#ifdef FLANN
            //#pragma omp critical (out)
#endif
            {

                // 3.4 loop over the candidates, do clipping
                for (set<int>::iterator it = neighborElems.begin(); it != neighborElems.end();
                        it++) {
                    int posSlaveNodes[N_P_E]; // the position of the node in the slave triangle
                    double slaveElemPrj[N_P_E * 3];

                    for (int ii = 0; ii < N_P_E; ii++) {
                        posSlaveNodes[ii] = slaveDirectElemTable[(*it) * N_P_E + ii];
                        for (int jj = 0; jj < 3; jj++) {
                            slaveElemPrj[ii * 3 + jj] = projections[posSlaveNodes[ii]][jj];
                        }
                    }
                    double result[N_P_E * N_P_E];
                    bool overlap = clip(masterElem, planeToProject, clipper, slaveElemPrj, result); // clip
                    if (overlap) { // only add the result to sparsity map if overlap happens, that is how C_BA is sparse
                        for (int ii = 0; ii < N_P_E; ii++) {
                            for (int jj = 0; jj < N_P_E; jj++) {
#ifdef FLANN
#pragma omp critical (map1)
#endif
                                {
                                    bool inserted = sparsityMapC_BA[posMasterNodes[ii]]->insert(
                                            pair<int, double>(posSlaveNodes[jj],
                                                    result[ii * N_P_E + jj])).second; // *.second is a bool indicating inserted or not
                                    if (!inserted)
                                        (sparsityMapC_BA[posMasterNodes[ii]]->at(posSlaveNodes[jj])) +=
                                                result[ii * N_P_E + jj]; // at() returns a reference, hence using "+=" is correct
                                } //omp critical
                            }
                        }
                    }
                    if (overlap && dual) { // if dual, compute sparsityMapC_BA_DUAL
                        double result_dual[N_P_E * N_P_E];
                        for (int ii = 0; ii < N_P_E * N_P_E; ii++)
                            result_dual[ii] = result[ii]; // copy the value in result
                        MortarMath::computeMatrixProduct(N_P_E, coeffMatrix, result_dual); // now it is dual
                        for (int ii = 0; ii < N_P_E; ii++) {
                            for (int jj = 0; jj < N_P_E; jj++) {
#ifdef FLANN
#pragma omp critical (map2)
#endif
                                {
                                    bool inserted =
                                            sparsityMapC_BA_DUAL[posMasterNodes[ii]]->insert(
                                                    pair<int, double>(posSlaveNodes[jj],
                                                            result_dual[ii * N_P_E + jj])).second; // *.second is a bool indicating inserted or not
                                    if (!inserted)
                                        (sparsityMapC_BA_DUAL[posMasterNodes[ii]]->at(
                                                posSlaveNodes[jj])) += result_dual[ii * N_P_E + jj]; // at() returns a reference, hence using "+=" is correct
                                } //omp critical
                            }
                        }
                    }
                }
            } //omp critical

            //delete
            for (map<int, double*>::iterator it = projections.begin(); it != projections.end();
                    it++)
                delete[] it->second;
        }

    } //#pragma omp parallel

    // 4. according to sparsity map, compute a, ia, ja of CSR formated C_BA
    C_BA_IA = new int[masterNumNodes + 1];
    C_BA_IA[0] = 1; // numbering from 1
    int nnz = 1; // number of non-zero entries
    for (int i = 0; i < masterNumNodes; i++) {
        nnz += sparsityMapC_BA[i]->size();
        C_BA_IA[i + 1] = nnz;
    }
    nnz--;
    C_BA_JA = new int[nnz];
    C_BA_A = new double[nnz];
    int count = 0;
    for (int i = 0; i < masterNumNodes; i++) {
        for (map<int, double>::iterator it = sparsityMapC_BA[i]->begin();
                it != sparsityMapC_BA[i]->end(); it++) {
            C_BA_JA[count] = it->first + 1; // numbering from 1
            C_BA_A[count] = it->second / factorC_BA;
            count++;
        }
    }assert(nnz == count);

    if (dual) {
        C_BA_A_DUAL = new double[nnz];
        int count = 0;
        for (int i = 0; i < masterNumNodes; i++) {
            for (map<int, double>::iterator it = sparsityMapC_BA_DUAL[i]->begin();
                    it != sparsityMapC_BA_DUAL[i]->end(); it++) {
                C_BA_A_DUAL[count] = it->second / factorC_BA;
                count++;
            }
        }
        assert(nnz == count);
    }

    // solve the problem that C_BB * 1 != C_BA * 1, this happens due to the projection of A not covering B
    if (!dual) {
        double *factor = new double[masterNumNodes];
        for (int i = 0; i < masterNumNodes; i++) {
            factor[i] = 0.0;
        }
        for (int i = 0; i < masterNumNodes; i++) {
            for (int j = C_BB_IA[i]; j < C_BB_IA[i + 1]; j++) {
                factor[i] += C_BB_A[j - 1];
            }
        }
        for (int i = 0; i < masterNumNodes; i++) {
            for (int j = C_BB_IA[i] + 1; j < C_BB_IA[i + 1]; j++) {
                factor[C_BB_JA[j - 1] - 1] += C_BB_A[j - 1];
            }
        }
        for (int i = 0; i < masterNumNodes; i++) {
            double sum = 0.0;
            for (int j = C_BA_IA[i]; j < C_BA_IA[i + 1]; j++) {
                sum += C_BA_A[j - 1];
            }assert(sum>factor[i]*1E-15);
            factor[i] /= sum;
        }
        for (int i = 0; i < masterNumNodes; i++) {
            for (int j = C_BA_IA[i]; j < C_BA_IA[i + 1]; j++) {
                C_BA_A[j - 1] *= factor[i];
            }
        }
        delete[] factor;
    } else {
        for (int i = 0; i < masterNumNodes; i++) {
            double factor = C_BB_A_DUAL[i];
            double sum = 0.0;
            for (int j = C_BA_IA[i]; j < C_BA_IA[i + 1]; j++) {
                sum += C_BA_A_DUAL[j - 1];
            }assert(sum>factor*1E-15);
            factor /= sum;
            for (int j = C_BA_IA[i]; j < C_BA_IA[i + 1]; j++) {
                C_BA_A_DUAL[j - 1] *= factor;
            }
        }
    }
    /*cout << endl << "=====C_BA_IA======" << endl;
     for (int i = 0; i < masterNumNodes + 1; i++)
     cout << C_BA_IA[i] << "   ";
     cout << endl << "=====C_BA_JA======" << endl;
     for (int i = 0; i < nnz; i++)
     cout << C_BA_JA[i] << "   ";
     cout << endl << "=====C_BA_A======" << endl;
     for (int i = 0; i < nnz; i++)
     cout << C_BA_A[i] << "   ";

     double aaa[4][4];
     int pp = 0;
     for (int i=0; i<4; i++)
     for (int j=0; j<4; j++) {
     aaa[i][j] = C_BA_A[pp];
     pp++;
     }

     double ttt[4] = {0,0,0,0};
     for (int i=0; i<4; i++)
     for (int j=0; j<4; j++)
     ttt[i]+=aaa[i][j];
     cout << endl << "++++++++ " << endl;
     for (int j=0; j<4; j++)
     cout << ttt[j] << endl;*/

    // delete
    for (int i = 0; i < masterNumNodes; i++)
        delete sparsityMapC_BA[i];
    delete[] sparsityMapC_BA;
    if (dual) {
        for (int i = 0; i < masterNumNodes; i++)
            delete sparsityMapC_BA_DUAL[i];
        delete[] sparsityMapC_BA_DUAL;
    }
}

void AbstractMortarMapper::initPardiso() {
#ifdef INTEL_MKL
    mtype = 2; // real symmetric matrix
    // set pardiso default parameters
    pardisoinit(pt, &mtype, iparm);
    iparm[2] = mklSetNumThreads;
    //cout << endl << "iparm[3]: " << iparm[2] << endl;
    maxfct = 1;// max number of factorizations
    mnum = 1;// which factorization to use
    msglvl = 0;// do NOT print statistical information
    neq = masterNumNodes;// number of rows of C_BB
    nrhs = 1;// number of right hand side
    int phase = 12;// analysis and factorization
    double ddum;// double dummy
    int idum;// integer dummy
    int error = 0;
    //cout<<"factorizing"<<endl;
    mkl_set_num_threads(mklSetNumThreads);
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, C_BB_A, C_BB_IA, C_BB_JA, &idum, &nrhs, iparm,
            &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        cerr << "Error in MortarMapper: pardiso factorization failed!" << error << endl;
        exit(EXIT_FAILURE);
    }
#endif
}

void AbstractMortarMapper::deletePardiso() {
#ifdef INTEL_MKL
    int phase = -1; // deallocate memory
    double ddum;// double dummy
    int idum;// integer dummy
    int error = 0;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &neq, C_BB_A, C_BB_IA, C_BB_JA, &idum, &nrhs, iparm,
            &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        cerr << "Error in MortarMapper: pardiso factorization failed!" << error << endl;
        exit(EXIT_FAILURE);
    }
#endif
}

void AbstractMortarMapper::initTables() {
    // using the map to store the nodeNumbers
    // but here the "key" is the node number, and the value is the position in nodeNumbers
    // the map is sorted automatically, so it is efficient for searching

    slaveDirectElemTable = new int[slaveNumElems * slaveNodesPerElem];
    masterDirectElemTable = new int[masterNumElems * masterNodesPerElem];
    slaveNodeToElemTable = new vector<int>*[slaveNumNodes];
    for (int i = 0; i < slaveNumNodes; i++)
        slaveNodeToElemTable[i] = new vector<int>;

    // 1. compute slaveDirectElemTable
    map<int, int> *slaveNodesMap = new map<int, int>();
    for (int i = 0; i < slaveNumNodes; i++)
        slaveNodesMap->insert(pair<int, int>(slaveNodeNumbers[i], i));
    for (int i = 0; i < slaveNumElems; i++)
        for (int j = 0; j < slaveNodesPerElem; j++)
            slaveDirectElemTable[i * slaveNodesPerElem + j] = slaveNodesMap->at(
                    slaveElemTable[i * slaveNodesPerElem + j]);
    delete slaveNodesMap;

    // 2. compute masterDirectElemTable
    map<int, int> *masterNodesMap = new map<int, int>();
    for (int i = 0; i < masterNumNodes; i++)
        masterNodesMap->insert(masterNodesMap->end(), pair<int, int>(masterNodeNumbers[i], i));
    for (int i = 0; i < masterNumElems; i++)
        for (int j = 0; j < masterNodesPerElem; j++)
            masterDirectElemTable[i * masterNodesPerElem + j] = masterNodesMap->at(
                    masterElemTable[i * masterNodesPerElem + j]);
    delete masterNodesMap;

    // 3. compute slaveNodeToElemTable
    for (int i = 0; i < slaveNumElems; i++)
        for (int j = 0; j < slaveNodesPerElem; j++) {
            int nodePos = slaveDirectElemTable[i * slaveNodesPerElem + j];
            slaveNodeToElemTable[nodePos]->push_back(i);
        }

    computeSlaveElemNormals();
}

void AbstractMortarMapper::deleteTables() {
    delete[] slaveDirectElemTable;
    delete[] masterDirectElemTable;
    for (int i = 0; i < slaveNumNodes; i++)
        delete slaveNodeToElemTable[i];
    delete[] slaveNodeToElemTable;
    delete[] slaveElemNormals;
}

void AbstractMortarMapper::initANNTree() {
#ifdef ANN
    ANNSlaveNodes = new double*[slaveNumNodes]; // ANN uses 2D array
    for (int i = 0; i < slaveNumNodes; i++)
    ANNSlaveNodes[i] = &slaveNodeCoors[i * 3];
    slaveNodesTree = new ANNkd_tree(ANNSlaveNodes, slaveNumNodes, 3);
#endif
#ifdef FLANN
    FLANNSlaveNodes = new flann::Matrix<double> (slaveNodeCoors, slaveNumNodes, 3);
    FLANNkd_tree = new flann::Index<flann::L2<double> >(*FLANNSlaveNodes, flann::KDTreeSingleIndexParams(1));
    FLANNkd_tree->buildIndex(); // Build binary tree for searching
#endif
}

void AbstractMortarMapper::deleteANNTree() {
#ifdef ANN
    delete[] ANNSlaveNodes;
    delete slaveNodesTree;
    annClose();
#endif
#ifdef FLANN
    delete FLANNSlaveNodes;
    delete FLANNkd_tree;
#endif
}

void AbstractMortarMapper::setN_P_E(int n) {
    assert(masterNodesPerElem == slaveNodesPerElem);
    assert(masterNodesPerElem == n);
    N_P_E = n;
}

bool AbstractMortarMapper::clip(const double *masterElem, int planeToProject,
        AbstractMortarMapper::PolygonPointClipper &clipper, const double *slaveElem,
        double *result) {
    // numbering the edge and point in the following way:
    // The edge i is the edge between point i and i+1.
    // The algorithm is from the book "J.Foley, Computer Graphics, Principles and Practice, 2nd edition, Addison-Wesley" P237 - P240.
    vector<double*> listInput;
    vector<double*> listOutput;
    double listInputStorage[20][3]; // 20 is enough to store end points of the clipped, increase it if the elements have more than 4 edges
    double listOutputStorage[20][3];
    for (int i = 0; i < N_P_E; i++) {
        copyPoint(&slaveElem[i * 3], listOutputStorage[i]);
        listOutput.push_back(listOutputStorage[i]);
    }
    //MortarMath::printElem(masterQuadPrj, 4);
    //MortarMath::printElem(slaveQuadPrj, 4);

    for (int i = 0; i < N_P_E; i++) {
        // i is the edge ID
        // 1. make the output the new input
        listInput.clear();
        for (unsigned j = 0; j < listOutput.size(); j++) {
            listInput.push_back(listInputStorage[listInput.size()]);
            double *newPoint = listInput.back();
            copyPoint(listOutput[j], newPoint);
        }
        listOutput.clear();
        const double *edgeP1 = &masterElem[i * 3];
        const double *edgeP2 = &masterElem[(i + 1) % N_P_E * 3];
        //MortarMath::printPoint(edgeP1);
        //MortarMath::printPoint(edgeP2);
        int size = listInput.size();

        // 2. clip the input list by edge, create new output list
        for (int j = 0; j < size; j++) {
            double *p1 = listInput[j];
            double *p2 = listInput[(j + 1) % size];
            bool inside1 = clipper.inside(i, p1);
            bool inside2 = clipper.inside(i, p2);
            if (!inside1 && inside2) {
                listOutput.push_back(listOutputStorage[listOutput.size()]);
                double *newPoint1 = listOutput.back();
                if (!MortarMath::intersect(edgeP1, edgeP2, p1, p2, planeToProject, newPoint1))
                    copyPoint(p1, newPoint1); // if two lines overlapped, make the out point the intersection
                listOutput.push_back(listOutputStorage[listOutput.size()]);
                double *newPoint2 = listOutput.back();
                copyPoint(p2, newPoint2);
            } else if (inside1 && inside2) {
                listOutput.push_back(listOutputStorage[listOutput.size()]);
                double *newPoint = listOutput.back();
                copyPoint(p2, newPoint);
            } else if (inside1 && !inside2) {
                listOutput.push_back(listOutputStorage[listOutput.size()]);
                double *newPoint = listOutput.back();
                if (!MortarMath::intersect(edgeP1, edgeP2, p1, p2, planeToProject, newPoint))
                    copyPoint(p2, newPoint); // if two lines overlapped, make the out point the intersection
            } else {
                // do nothing
            }
        }
    }

    // remove the overlapped points in the polygon
    double nonOverlapStorage[N_P_E * 2][3];
    const double EPS = 1E-15;
    const double TOL_SQR = EPS * EPS * longestEdgeLengthSquare(masterElem, N_P_E);
    int count = 0;
    int size = listOutput.size();
    for (int j = 0; j < size; j++) {
        double *p1 = listOutput[j];
        double *p2 = listOutput[(j + 1) % size];
        if (MortarMath::distanceSquare(p1, p2) > TOL_SQR) {
            copyPoint(p2, nonOverlapStorage[count]);
            count++;
        }
    }
    listOutput.clear();
    for (int j = 0; j < count; j++)
        listOutput.push_back(nonOverlapStorage[j]);

    if (listOutput.size() < 3)
        return false;

    gaussQuadratureOnClip(masterElem, slaveElem, planeToProject, listOutput, result);
    return true;
}

void AbstractMortarMapper::gaussQuadratureOnClip(const double *elem1, const double *elem2,
        int planeToProject, vector<double*> &polygon, double *result) {
    for (int i = 0; i < N_P_E * N_P_E; i++)
        result[i] = 0;

    // if the clipped polygon is a triangle, do Gauss quadrature on it directly
    if (polygon.size() == 3) {
        double clipTriangle[9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                clipTriangle[i * 3 + j] = polygon[i][j];

        Integrand *integrand = getShapeFuncProductOn2Elems(elem1, elem2, planeToProject,
                clipTriangle);

        for (int i = 0; i < N_P_E; i++) {
            for (int j = 0; j < N_P_E; j++) {
                integrand->setShapeFunctions(i, j);
                result[i * N_P_E + j] = GaussQuadrature::gaussQuadratureOnTriangle(integrand);
            }
        }
        delete integrand;
    } else { // if not, divide it into triangles, on each of which doing Gauss quadrature
        double center[3];
        polygonCenter(polygon, center);
        int size = polygon.size();
        for (int i = 0; i < size; i++) {
            double clipTriangle[9];
            buildTrianagle(center, polygon[i], polygon[(i + 1) % size], clipTriangle);
            Integrand *integrand = getShapeFuncProductOn2Elems(elem1, elem2, planeToProject,
                    clipTriangle);
            for (int i = 0; i < N_P_E; i++) {
                for (int j = 0; j < N_P_E; j++) {
                    integrand->setShapeFunctions(i, j);
                    result[i * N_P_E + j] += GaussQuadrature::gaussQuadratureOnTriangle(integrand);
                }
            }
            delete integrand;
        }
    }
}

double AbstractMortarMapper::computeSearchRadiusSquare(const double* masterElem) {
    // 1. find the longest edge of this element (masterElem)
    double masterElemCopy[3 * N_P_E];
    copyElem(masterElem, N_P_E, masterElemCopy);
    double lengthSqr = longestEdgeLengthSquare(masterElemCopy, N_P_E);

    // 2. find the longest edge of the neighboring elements (the slave side)
    // the neighboring elements are found by
    // first find all neighboring nodes of the nodes of the master element,
    // then all elements on the slave side having these nodes are the neighboring elements.
    // These elements together with the master element itself help to compute a reasonable search radius.
    double dummy;
    int nb;
    vector<int>* elemVec[N_P_E];
    for (int i = 0; i < N_P_E; i++) {
#ifdef ANN
        slaveNodesTree->annkSearch(&masterElemCopy[i * 3], 1, &nb, &dummy);
#endif
#ifdef FLANN
        flann::Matrix<double> masterElemCopyFlann(&masterElemCopy[i*3], 1, 3);
        vector<vector<int> > indexes_tmp;
        vector<vector<double> > dists_tmp;
        FLANNkd_tree->knnSearch(masterElemCopyFlann, indexes_tmp, dists_tmp, 1, flann::SearchParams(1));
        nb = indexes_tmp[0].size();
#endif

        elemVec[i] = slaveNodeToElemTable[nb]; // do not modify/delete contents in elemVec!!!
    }

    double searchRadiusSqr = lengthSqr; // the goal of the first two steps is to set up this value
    set<int> neighborElemsTmp; // set allows no multiple entries
    for (int i = 0; i < N_P_E; i++)
        for (unsigned j = 0; j < elemVec[i]->size(); j++)
            neighborElemsTmp.insert(elemVec[i]->at(j));

    double elemTmp[3 * N_P_E];
    for (set<int>::iterator it = neighborElemsTmp.begin(); it != neighborElemsTmp.end(); it++) {
        // for a single element, find its longest edge
        getElemCoor(*it, AbstractMortarMapper::SLAVE, elemTmp);
        lengthSqr = longestEdgeLengthSquare(elemTmp, N_P_E);
        if (lengthSqr > searchRadiusSqr)
            searchRadiusSqr = lengthSqr;
    }

    // 3. the radius is 1.12 times the longest edge in the first two steps
    searchRadiusSqr *= 1.12;
    searchRadiusSqr *= 1.12;
    return searchRadiusSqr;
}

void AbstractMortarMapper::findCandidates(const double* masterElem, const double* masterElemNormal,
        double searchRadiusSqr, set<int> &neighborElems) {
    // Use fixed radius search on end points of the element,
    // all elements containing these points are the overlapped candidates
    // 1. find all neighboring elements in the radius
    double masterElemCopy[N_P_E * 3];
    copyElem(masterElem, N_P_E, masterElemCopy);
    // OpenMP parallelize this loop
    // Ann is not thread safe
    for (int i = 0; i < N_P_E; i++) {
#ifdef ANN
        int n_nb = 0;
        n_nb = slaveNodesTree->annkFRSearch(&masterElemCopy[i * 3], searchRadiusSqr, 0); // get the number of neighbors in a radius
        int nbs[n_nb];
        double dummy[n_nb];
        slaveNodesTree->annkFRSearch(&masterElemCopy[i * 3], searchRadiusSqr, n_nb, nbs, dummy);// get the real neighbors (ANN uses the square of the radius)
#endif

#ifdef FLANN
        flann::Matrix<double> masterElemCopyFlann(&masterElemCopy[i*3], 1, 3);
        vector<vector<int> > indexes_tmp;
        vector<vector<double> > dists_tmp;
        FLANNkd_tree->radiusSearch(masterElemCopyFlann, indexes_tmp, dists_tmp, searchRadiusSqr, flann::SearchParams(1));
        int n_nb = indexes_tmp[0].size();
        int* nbs = &indexes_tmp[0][0];
#endif
        for (int j = 0; j < n_nb; j++) {
            vector<int> * vecTmp = slaveNodeToElemTable[nbs[j]];
            for (unsigned k = 0; k < vecTmp->size(); k++) {
                neighborElems.insert(vecTmp->at(k));
            }
        }
    }
    vector<set<int>::iterator> toDelete;
    // 2. kick out the elements that have wrong normal direction
    for (set<int>::iterator it = neighborElems.begin(); it != neighborElems.end(); it++) {
        double *slaveElemNormal = &slaveElemNormals[(*it) * 3];
        if (kickOutCandidate(masterElemNormal, slaveElemNormal, 0.7))
            toDelete.push_back(it); //neighborElems.erase(it) is wrong! because it++ goes to a magic place later
    }
    if (toDelete.size() == neighborElems.size()) { // loose the restriction if no neighbors are found
        toDelete.clear();
        for (set<int>::iterator it = neighborElems.begin(); it != neighborElems.end(); it++) {
            double *slaveElemNormal = &slaveElemNormals[(*it) * 3];
            if (kickOutCandidate(masterElemNormal, slaveElemNormal, 0.01))
                toDelete.push_back(it); //neighborElems.erase(it) is wrong! because it++ goes to a magic place later
        }
    }
    for (int i = 0; i < toDelete.size(); i++)
        neighborElems.erase(toDelete[i]); // out kicking
    assert(neighborElems.size()!=0);
}

bool AbstractMortarMapper::kickOutCandidate(const double *masterUnitNormal,
        const double *slaveUnitNormal, double bound) {
    //acos(0.7) = pi/4
    //acos(0.01) = pi/2
    if (!oppositeSurfaceNormal)
        return (MortarMath::dotProduct(masterUnitNormal, slaveUnitNormal) < bound);
    return (MortarMath::dotProduct(masterUnitNormal, slaveUnitNormal) > -bound);
}

void AbstractMortarMapper::projectToElemPlane(const double *elem, const double *planeUnitNormal,
        set<int> &neighborElems, map<int, double*> &projections) {
    set<int> neighborNodes;
    for (set<int>::iterator it = neighborElems.begin(); it != neighborElems.end(); it++) {
        for (int i = 0; i < N_P_E; i++)
            neighborNodes.insert(slaveDirectElemTable[(*it) * N_P_E + i]);
    }
    //cout << neighborNodes.size();
    for (set<int>::iterator it = neighborNodes.begin(); it != neighborNodes.end(); it++) {
        double *projectionTmp = new double[3];
        MortarMath::projectToElemPlane(&elem[0], planeUnitNormal, &slaveNodeCoors[(*it) * 3], 1,
                projectionTmp);
        projections.insert(projections.end(), pair<int, double*>(*it, projectionTmp));
    }
}

void AbstractMortarMapper::getElemCoor(int elemIndex, MeshLabel label, double *elem) {
    // compute the coordinates of a triangle by its id
    int *elemTable;
    int n_p_e;
    double *nodeCoors;
    if (label == MASTER) {
        elemTable = masterDirectElemTable;
        n_p_e = masterNodesPerElem;
        nodeCoors = masterNodeCoors;
    } else {
        elemTable = slaveDirectElemTable;
        n_p_e = slaveNodesPerElem;
        nodeCoors = slaveNodeCoors;
    }

    for (int i = 0; i < n_p_e; i++) {
        int nodePos = elemTable[elemIndex * n_p_e + i]; // position of the node
        for (int j = 0; j < 3; j++)
            elem[i * 3 + j] = nodeCoors[nodePos * 3 + j];
    }
}

void AbstractMortarMapper::copyPoint(const double *origin, double *copy) {
    for (int i = 0; i < 3; i++)
        copy[i] = origin[i];
}

void AbstractMortarMapper::copyElem(const double *origin, int size, double *copy) {
    for (int i = 0; i < 3 * size; i++)
        copy[i] = origin[i];
}

void AbstractMortarMapper::buildTrianagle(const double *p0, const double *p1, const double *p2,
        double *triangle) {
    for (int i = 0; i < 3; i++) {
        triangle[i] = p0[i];
        triangle[3 + i] = p1[i];
        triangle[6 + i] = p2[i];
    }
}

void AbstractMortarMapper::polygonCenter(vector<double*> &polygon, double *center) {
    for (int j = 0; j < 3; j++)
        center[j] = 0;
    for (unsigned i = 0; i < polygon.size(); i++)
        for (int j = 0; j < 3; j++)
            center[j] += polygon[i][j];

    int size = polygon.size();
    for (int j = 0; j < 3; j++)
        center[j] /= size;
}

double AbstractMortarMapper::longestEdgeLengthSquare(const double *elem, int size) {
    double l[size];
    for (int i = 0; i < size; i++) {
        l[i] = MortarMath::distanceSquare(&elem[i * 3], &elem[(i + 1) % size * 3]);
    }
    double max = l[0];
    for (int i = 1; i < size; i++) {
        if (l[i] > max)
            max = l[i];
    }
    return max;
}

void AbstractMortarMapper::computeDualCoeffMatrix(const double *elem, double *coeffMatrix) {
    double massMatrixLower[N_P_E * N_P_E]; // lower part of the symmetric matrix
    Integrand *integrand = getShapeFuncProductOnSingleElem(elem);
    for (int i = 0; i < N_P_E; i++) {
        for (int j = i; j < N_P_E; j++) {
            integrand->setShapeFunctions(i, j);
            if (N_P_E == 4)
                massMatrixLower[i * N_P_E + j] = GaussQuadrature::gaussQuadratureOnQuad(integrand);
            else if (N_P_E == 3)
                massMatrixLower[i * N_P_E + j] = GaussQuadrature::gaussQuadratureOnTriangle(
                        integrand);
        }
    }
    double shapeFuncIntegration[N_P_E];
    for (int i = 0; i < N_P_E; i++) {
        integrand->setShapeFunctions(i, -1);
        if (N_P_E == 4)
            shapeFuncIntegration[i] = GaussQuadrature::gaussQuadratureOnQuad(integrand);
        else if (N_P_E == 3)
            shapeFuncIntegration[i] = GaussQuadrature::gaussQuadratureOnTriangle(integrand);
    }
    for (int i = 0; i < N_P_E; i++) {
        for (int j = 0; j < N_P_E; j++) {
            coeffMatrix[i * N_P_E + j] = 0.0;
        }
    }
    for (int i = 0; i < N_P_E; i++)
        coeffMatrix[i * N_P_E + i] = shapeFuncIntegration[i];
    delete integrand;
    int dummy[N_P_E];
    int info = -1;
#ifdef LAPACK
#ifdef INTEL_MKL
    mkl_set_num_threads(mklSetNumThreads);
#endif
    info = LAPACKE_dsysv(LAPACK_COL_MAJOR, 'L', N_P_E, N_P_E, massMatrixLower, N_P_E, dummy,
            coeffMatrix, N_P_E);
#endif
    assert(info == 0);
}

void AbstractMortarMapper::checkNullPointers() {
    if (!dual) {
        assert(C_BB_A!=NULL);
        assert(C_BB_IA!=NULL);
        assert(C_BB_JA!=NULL);
        assert(C_BB_A_DUAL==NULL);

        assert(C_BA_A!=NULL);
        assert(C_BA_IA!=NULL);
        assert(C_BA_JA!=NULL);
        assert(C_BA_A_DUAL==NULL);
    } else {
        assert(C_BB_A==NULL);
        assert(C_BB_IA==NULL);
        assert(C_BB_JA==NULL);
        assert(C_BB_A_DUAL!=NULL);
        // if it is dual, we still have to store C_BA in order to do press2force mapping
        assert(C_BA_A!=NULL);
        assert(C_BA_IA!=NULL);
        assert(C_BA_JA!=NULL);
        assert(C_BA_A_DUAL!=NULL);
    }
}
