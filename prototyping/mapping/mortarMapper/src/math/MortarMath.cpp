#include "MortarMath.h"
#include <math.h>
#include <mkl_lapacke.h>
#include <iostream>
#include <stdlib.h>
#include <assert.h>

using namespace std;
double* GaussQuadrature::triGaussPoints3;
double* GaussQuadrature::triWeights3;
double* GaussQuadrature::triGaussPoints6;
double* GaussQuadrature::triWeights6;
double* GaussQuadrature::triGaussPoints7;
double* GaussQuadrature::triWeights7;
double* GaussQuadrature::triGaussPoints12;
double* GaussQuadrature::triWeights12;
double* GaussQuadrature::quadGaussPoints1;
double* GaussQuadrature::quadWeights1;
double* GaussQuadrature::quadGaussPoints4;
double* GaussQuadrature::quadWeights4;
double* GaussQuadrature::quadGaussPoints9;
double* GaussQuadrature::quadWeights9;

double GaussQuadrature::gaussQuadratureOnTriangle(
		Integrand *integrand) {

	int numGaussPoints = integrand->getNumGaussPoints();
	const double *weights = getTriangleWeights(numGaussPoints);

	double toReturn = 0;
	for (int i = 0; i < numGaussPoints; i++)
		toReturn += weights[i] * integrand->evaluate(i);
	return toReturn;
}

double GaussQuadrature::gaussQuadratureOnQuad(
		Integrand *integrand) {
	int numGaussPoints = integrand->getNumGaussPoints();
	const double *weights = getQuadWeights(numGaussPoints);

	double toReturn = 0;
	for (int i = 0; i < numGaussPoints; i++)
		toReturn += weights[i] * integrand->evaluate(i);
	return toReturn;
}

const double* GaussQuadrature::getTriangleGaussPoints(
		int numGaussPoints) {
	if (numGaussPoints == 3) {
		if (triGaussPoints3 == NULL) {
			triGaussPoints3 = new double[9];
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					if (i==j)
						triGaussPoints3[i*3+j] = 2.0 / 3.0;
					else
						triGaussPoints3[i*3+j] = 1.0 / 6.0;
		}
		return triGaussPoints3;
	} else if (numGaussPoints == 6) {
		if (triGaussPoints6 == NULL) {
			triGaussPoints6 = new double[18];
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					if (i==j)
						triGaussPoints6[i*3+j] = 0.816847572980459;
					else
						triGaussPoints6[i*3+j] = 0.091576213509771;
			for (int i=3; i<6; i++)
				for (int j=0; j<3; j++)
					if (i-3==j)
						triGaussPoints6[i*3+j] = 0.108103018168070;
					else
						triGaussPoints6[i*3+j] = 0.445948490915965;
		}
		return triGaussPoints6;
	} else if (numGaussPoints == 7) {
		if (triGaussPoints7 == NULL) {
			triGaussPoints7 = new double[21];
			for (int i=0; i<3; i++)
				triGaussPoints7[i] = 1.0 / 3.0;
			for (int i=1; i<4; i++)
				for (int j=0; j<3; j++)
					if (i-1==j)
						triGaussPoints7[i*3+j] = 0.797426985353087;
					else
						triGaussPoints7[i*3+j] = 0.101286507323456;
			for (int i=4; i<7; i++)
				for (int j=0; j<3; j++)
					if (i-4==j)
						triGaussPoints7[i*3+j] = 0.059715871789770;
					else
						triGaussPoints7[i*3+j] = 0.470142064105115;
		}
		return triGaussPoints7;
	} else if (numGaussPoints == 12) {
		if (triGaussPoints12 == NULL) {
			triGaussPoints12 = new double[36];
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					if (i==j)
						triGaussPoints12[i*3+j] = 0.873821971016996;
					else
						triGaussPoints12[i*3+j] = 0.063089014491502;
			for (int i=3; i<6; i++)
				for (int j=0; j<3; j++)
					if (i-3==j)
						triGaussPoints12[i*3+j] = 0.501426509658179;
					else
						triGaussPoints12[i*3+j] = 0.249286745170910;

			double v[3] = {0.636502499121399, 0.310352451033785, 0.053145049844816};
			for (int i=0; i<3; i++) {
				triGaussPoints12[6*3+i] = v[i];
				triGaussPoints12[7*3+i] = v[(i+1) % 3];
				triGaussPoints12[8*3+i] = v[(i+2) % 3];
			}
			for (int i=0; i<3; i++) {
				triGaussPoints12[9*3+i] = v[(-i+3)%3];
				triGaussPoints12[10*3+i] = v[(-i+1+3) % 3];
				triGaussPoints12[11*3+i] = v[(-i+2+3) % 3];
			}
		}
		return triGaussPoints12;
	} else {
		cerr << "Error in MortarMath::GaussQuadrature::getTriangleGaussPoints()" << endl;
		cerr << "\t wrong number of Gauss points" << endl;
		exit(EXIT_FAILURE);
	}
}

const double* GaussQuadrature::getTriangleWeights(
		int numGaussPoints) {
	if (numGaussPoints == 3) {
		if (triWeights3 == NULL) {
			triWeights3 = new double[3];
			for (int i=0; i<3; i++)
				triWeights3[i] = 1.0 / 3.0;
		}
		return triWeights3;
	} else if (numGaussPoints == 6) {
		if (triWeights6 == NULL) {
			triWeights6 = new double[6];
			for (int i=0; i<3; i++)
				triWeights6[i] = 0.109951743655322;
			for (int i=3; i<6; i++)
				triWeights6[i] = 0.223381589678011;
		}
		return triWeights6;
	} else if(numGaussPoints == 7) {
		if (triWeights7 == NULL) {
			triWeights7 = new double[7];
			triWeights7[0] = 0.225000000000000;
			for (int i=1; i<4; i++)
				triWeights7[i] = 0.125939180544827;
			for (int i=4; i<7; i++)
				triWeights7[i] = 0.132394152788506;
		}
		return triWeights7;
	} else if(numGaussPoints == 12) {
		if (triWeights12 == NULL) {
			triWeights12 = new double[12];
			for (int i=0; i<3; i++)
				triWeights12[i] = 0.050844906370207;
			for (int i=3; i<6; i++)
				triWeights12[i] = 0.116786275726379;
			for (int i=6; i<12; i++)
				triWeights12[i] = 0.082851075618374;
		}
		return triWeights12;
	} else {
		cerr << "Error in MortarMath::GaussQuadrature::getTriangleWeights()" << endl;
		cerr << "\t wrong number of Gauss points" << endl;
		exit(EXIT_FAILURE);
	}
}

const double* GaussQuadrature::getQuadGaussPoints(
		int numGaussPoints) {
	if (numGaussPoints == 1) {
		if (quadGaussPoints1 == NULL) {
			quadGaussPoints1 = new double[2];
			quadGaussPoints1[0] = 0.0;
			quadGaussPoints1[1] = 0.0;
		}
		return quadGaussPoints1;
	} else if (numGaussPoints == 4) {
		if (quadGaussPoints4 == NULL) {
			quadGaussPoints4 = new double[8];
			double v = sqrt(3.0) / 3.0;
			quadGaussPoints4[0] = v;
			quadGaussPoints4[1] = v;
			quadGaussPoints4[2] = -v;
			quadGaussPoints4[3] = v;
			quadGaussPoints4[4] = -v;
			quadGaussPoints4[5] = -v;
			quadGaussPoints4[6] = v;
			quadGaussPoints4[7] = -v;
		}
		return quadGaussPoints4;
	} else if (numGaussPoints == 9) {
		if (quadGaussPoints9 == NULL) {
			quadGaussPoints9 = new double[18];
			double v = sqrt(3.0 / 5.0);
			quadGaussPoints9[0] = 0;
			quadGaussPoints9[1] = 0;

			quadGaussPoints9[2] = -v;
			quadGaussPoints9[3] = 0;
			quadGaussPoints9[4] = 0;
			quadGaussPoints9[5] = -v;
			quadGaussPoints9[6] = v;
			quadGaussPoints9[7] = 0;
			quadGaussPoints9[8] = 0;
			quadGaussPoints9[9] = v;

			quadGaussPoints9[10] = -v;
			quadGaussPoints9[11] = -v;
			quadGaussPoints9[12] = v;
			quadGaussPoints9[13] = -v;
			quadGaussPoints9[14] = -v;
			quadGaussPoints9[15] = v;
			quadGaussPoints9[16] = v;
			quadGaussPoints9[17] = v;
		}
		return quadGaussPoints9;
	} else {
		cerr << "Error in MortarMath::GaussQuadrature::getQuadGaussPoints()" << endl;
		cerr << "\t wrong number of Gauss points" << endl;
		exit(EXIT_FAILURE);
	}
}

const double* GaussQuadrature::getQuadWeights(
		int numGaussPoints) {
	if (numGaussPoints == 1) {
		if (quadWeights1 == NULL) {
			quadWeights1 = new double[1];
			quadWeights1[0] = 4.0;
		}
		return quadWeights1;
	} else if (numGaussPoints == 4) {
		if (quadWeights4 == NULL) {
			quadWeights4 = new double[4];
			for (int i=0; i<4; i++)
				quadWeights4[i] = 1.0;
		}
		return quadWeights4;
	} else if (numGaussPoints == 9) {
		if (quadWeights9 == NULL) {
			quadWeights9 = new double[9];
			quadWeights9[0] = 8.0/9.0 * 8.0/9.0;
			for (int i=1; i<5; i++)
				quadWeights9[i] = 8.0/9.0 * 5.0/9.0;
			for (int i=5; i<9; i++)
				quadWeights9[i] = 5.0/9.0 * 5.0/9.0;
		}
		return quadWeights9;
	} else {
		cerr << "Error in MortarMath::GaussQuadrature::getQuadWeights()" << endl;
		cerr << "\t wrong number of Gauss points" << endl;
		exit(EXIT_FAILURE);
	}
}

void MortarMath::projectToElemPlane(const double* pointOnPlane,
		const double *unitNormal, const double *points,
		int num, double *projections) {

	double ridge[3];
	double distance;
	double path[3];
	for (int i=0; i<num; i++) {
		for (int j=0; j<3; j++)
			ridge[j] = pointOnPlane[j] - points[i*3+j];
		distance = dotProduct(unitNormal, ridge);
		for (int j=0; j<3; j++)
			path[j] = distance * unitNormal[j];

		for (int j=0; j<3; j++)
			projections[i*3 + j] = points[i*3+j] + path[j];
	}
}

void MortarMath::computeTriangleNormal(const double *triangle, bool unitLength,
		double *normal) {
	const int NUM = 3; // three points in a triangle
	double length = 0; // length of the normal vector
	normal[0] = (triangle[2 * NUM + 1] - triangle[0 * NUM + 1]) * (triangle[2
			* NUM + 2] - triangle[1 * NUM + 2]) - (triangle[2 * NUM + 1]
			- triangle[1 * NUM + 1]) * (triangle[2 * NUM + 2] - triangle[0
			* NUM + 2]);
	normal[1] = -((triangle[2 * NUM + 0] - triangle[0 * NUM + 0]) * (triangle[2
			* NUM + 2] - triangle[1 * NUM + 2]) - (triangle[2 * NUM + 0]
			- triangle[1 * NUM + 0]) * (triangle[2 * NUM + 2] - triangle[0
			* NUM + 2]));
	normal[2] = (triangle[2 * NUM + 0] - triangle[0 * NUM + 0]) * (triangle[2
			* NUM + 1] - triangle[1 * NUM + 1]) - (triangle[2 * NUM + 0]
			- triangle[1 * NUM + 0]) * (triangle[2 * NUM + 1] - triangle[0
			* NUM + 1]);

	if (!unitLength)
		return;

	length = vectorLength(normal);

	for (int i = 0; i < 3; i++)
		normal[i] /= length;
}

double MortarMath::computeArea(const double *triangle) {
	double area = 0.0;
	double normal[3];
	computeTriangleNormal(triangle, false, normal); // 2 times the length of the normal vector equals the area
	area = vectorLength(normal);
	area = 0.5 * area;
	return area;
}

bool MortarMath::computeLocalCoorInTriangle(const double *triangle, int planeToProject, const double *point,
		double *localCoor) {
	/* Normally, the local coordinates can be solved by:
	 *  |x_1 x_2 x_3|   |xi_1|    |x_4|
	 *  |y_1 y_2 y_3| * |xi_2| =  |y_4|
	 *  |z_1 z_2 z_3|   |xi_3|    |z_4|
	 *
	 * But if x_1 = x_2 = x_3 = 0 than the system cannot be solved.
	 * So we remove one equation by xi_1 + xi_2 + xi_3 = 1.
	 * This indicates projection.
	 * Choose among planes (xy or yz or zx) the one which has the smallest angle
	 * with the triangle normal.
	 * For example, if the angle between of normals of xy-plane and the triangle
	 * is 0, then the system is replaced by
	 *  |1.0 1.0 1.0|   |xi_1|    |1.0|
	 *  |y_1 y_2 y_3| * |xi_2| =  |y_4|
	 *  |z_1 z_2 z_3|   |xi_3|    |z_4|
	 */
	double A[9];
	for (int i = 0; i < 9; i++)
		A[i] = triangle[i];

	for (int i = 0; i < 3; i++)
		localCoor[i] = point[i];

	A[planeToProject] = A[planeToProject + 3] = A[planeToProject + 6]
	                                              = localCoor[planeToProject] = 1.0;

	/*int dummy[3];
	lapack_int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 3, 1, A, 3, dummy,
			localCoor, 3);
	if (info != 0) {
		cerr << "ERROR:\t triangle interpolation failed!" << endl;
		exit(EXIT_FAILURE);
	}*/
	solve3x3LinearSystem(A, planeToProject, localCoor);
	for (int i=0; i<3; i++) {
		if (localCoor[i] > 1.0)
			return false;
		if (localCoor[i] < 0)
			return false;
	}
	// make sure the sum is 1.0
	//assert(fabs(localCoor[0] + localCoor[1] + localCoor[2] -1) < 1E-15);
	localCoor[0] = 1.0 - localCoor[1] - localCoor[2];
	return true;
}

void MortarMath::computeGlobalCoorInTriangle(const double *triangle,
		const double *localCoor, double *globalCoor) {
	for (int i = 0; i < 3; i++) {
		globalCoor[i] = 0;
		for (int j = 0; j < 3; j++)
			globalCoor[i] += localCoor[j] * triangle[i + j * 3];
	}
}

double MortarMath::getQuadShapeFunc(int id, const double *xi_eta) {
	if (id == 0) {
		return 0.25 * (1.0 - xi_eta[0]) * (1.0 - xi_eta[1]);
	} else if (id == 1) {
		return 0.25 * (1.0 + xi_eta[0]) * (1.0 - xi_eta[1]);
	} else if (id == 2) {
		return 0.25 * (1.0 + xi_eta[0]) * (1.0 + xi_eta[1]);
	} else if (id == 3) {
		return 0.25 * (1.0 - xi_eta[0]) * (1.0 + xi_eta[1]);
	}
	return -100000000.0;
}

double MortarMath::computeQuadDetJ(const double* quad, const double *xi_eta) {
	// d_N_d_xi[4] contains the partial derivative w.r.t. xi of the shape functions
	double d_N_d_xi[4];
	// d_N_d_eta[4] contains the partial derivative w.r.t. eta of the shape functions
	double d_N_d_eta[4];

	d_N_d_xi[0] = - 0.25 * (1 - xi_eta[1]);
	d_N_d_xi[1] = - d_N_d_xi[0];
	d_N_d_xi[2] = 0.25 * (1 + xi_eta[1]);
	d_N_d_xi[3] = - d_N_d_xi[2];

	d_N_d_eta[0] = - 0.25 * (1 - xi_eta[0]);
	d_N_d_eta[1] = - 0.25 * (1 + xi_eta[0]);
	d_N_d_eta[2] = - d_N_d_eta[1];
	d_N_d_eta[3] = - d_N_d_eta[0];

	// g1 and g2 are the local basis vectors, and det(J)=||g1 x g2||
	double g1[3] = {0,0,0};
	double g2[3] = {0,0,0};

	for (int i=0; i<3; i++) {
		for (int j=0; j<4; j++) {
			g1[i] += quad[j*3 + i] * d_N_d_xi[j];
			g2[i] += quad[j*3 + i] * d_N_d_eta[j];
		}
	}

	double crossProduct[3];
	crossProduct[0] = g1[1]*g2[2] - g2[1]*g1[2];
	crossProduct[1] = -(g1[0]*g2[2] - g2[0]*g1[2]);
	crossProduct[2] = g1[0]*g2[1] - g2[0]*g1[1];
	return vectorLength(crossProduct);
}

void MortarMath::computeQuadNormal(const double *quad, bool unitLength,
		double *normal) {
	// d_N_d_xi[4] contains the partial derivative w.r.t. xi of the shape functions at (0, 0)
	double d_N_d_xi[4];
	// d_N_d_eta[4] contains the partial derivative w.r.t. eta of the shape functions at (0, 0)
	double d_N_d_eta[4];

	d_N_d_xi[0] = - 0.25;
	d_N_d_xi[1] = 0.25;
	d_N_d_xi[2] = 0.25;
	d_N_d_xi[3] = - 0.25;

	d_N_d_eta[0] = - 0.25;
	d_N_d_eta[1] = - 0.25;
	d_N_d_eta[2] = 0.25;
	d_N_d_eta[3] = 0.25;
	// g1 and g2 are the local basis vectors, and det(J)=||g1 x g2||
	double g1[3] = {0,0,0};
	double g2[3] = {0,0,0};

	for (int i=0; i<3; i++) {
		for (int j=0; j<4; j++) {
			g1[i] += quad[j*3 + i] * d_N_d_xi[j];
			g2[i] += quad[j*3 + i] * d_N_d_eta[j];
		}
	}
	double length;
	normal[0] = g1[1]*g2[2] - g2[1]*g1[2];
	normal[1] = -(g1[0]*g2[2] - g2[0]*g1[2]);
	normal[2] = g1[0]*g2[1] - g2[0]*g1[1];

	if (!unitLength)
		return;

	length = vectorLength(normal);

	for (int i = 0; i < 3; i++)
		normal[i] /= length;
}

void MortarMath::getQuadCenter(const double *quad, double *center) {
	for (int i=0; i<3; i++)
		center[i] = 0;
	for (int i=0; i<4; i++)
		for (int j=0; j<3; j++)
			center[j] += quad[i*3 + j];
	for (int i=0; i<3; i++)
		center[i] /= 4.0;
}

bool MortarMath::computeLocalCoorInQuad(const double *quad, int planeToProject,
		const double *point, double *localCoor) {
	/*
	 * So we use two coordinates among x, y, z.
	 * This indicates projection.
	 * Choose among planes (xy or yz or zx) the one which has the smallest angle
	 * with the quad normal.
	 */
	double x[4];
	double y[4];
	int x_direc = (planeToProject + 1) % 3;
	int y_direc = (planeToProject + 2) % 3;
	for (int i=0; i<4; i++) {
		x[i] = quad[i*3 + x_direc];
		y[i] = quad[i*3 + y_direc];
	}
	double x0 = point[x_direc];
	double y0 = point[y_direc];

	double a1 = x[0]+x[1]+x[2]+x[3] - 4.0*x0;
	double b1 = -x[0]+x[1]+x[2]-x[3];
	double c1 = -x[0]-x[1]+x[2]+x[3];
	double d1 = x[0]-x[1]+x[2]-x[3];

	double a2 = y[0]+y[1]+y[2]+y[3] - 4.0*y0;
	double b2 = -y[0]+y[1]+y[2]-y[3];
	double c2 = -y[0]-y[1]+y[2]+y[3];
	double d2 = y[0]-y[1]+y[2]-y[3];

	double delta[2];
	double J_T[4]; // transpose of Jacobian --- to use column major in lapack
	double F[2]; // -F
	localCoor[0] = 0;
	localCoor[1] = 0;
	//int dummy[2];
	const double EPS = 1E-15;

	const int MAX_ITER_NUM = 100;
	for (int i=0; i<MAX_ITER_NUM; i++) {
		J_T[0] = b1 + d1 * localCoor[1];
		J_T[2] = c1 + d1 * localCoor[0];
		J_T[1] = b2 + d2 * localCoor[1];
		J_T[3] = c2 + d2 * localCoor[0];
		F[0] = a1 + b1 * localCoor[0] + c1 * localCoor[1] + d1 * localCoor[0] * localCoor[1];
		F[1] = a2 + b2 * localCoor[0] + c2 * localCoor[1] + d2 * localCoor[0] * localCoor[1];
		delta[0] = -F[0];
		delta[1] = -F[1];
		/*int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 2, 1, J_T, 2, dummy,
				delta, 2);
		if (info != 0) {
			cerr << "ERROR in MortarMath::computeLocalCoorInQuad()!" << endl;
			exit(EXIT_FAILURE);
		}*/
		solve2x2LinearSystem(J_T, delta);
		if (fabs(delta[0]) < EPS && fabs(delta[1]) < EPS) {
			assert(i < 10);
			break;
		}
		localCoor[0] += delta[0];
		localCoor[1] += delta[1];
	}
	for (int i=0; i<2; i++) {
		if (localCoor[i] > 1.0)
			return false;
		if (localCoor[i] < -1.0)
			return false;
	}
	return true;
}

double MortarMath::distanceSquare(const double *p1, const double *p2) {
	double d;
	d = (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1])
			+ (p1[2] - p2[2]) * (p1[2] - p2[2]);
	return d;
}

bool MortarMath::intersect(const double *la0, const double *la1,
		const double *lb0, const double *lb1, int planeToProject,
		double *intersection) {

	// if lb0 and lb1 are overlapped, big numerical error for intersection happens
	const double EPS = 1E-6;
	const double TOL_SQR = EPS * EPS * distanceSquare(la0, la1);
	if (distanceSquare(lb0, lb1) < TOL_SQR)
		return false;
	// equation (in vector): alpha a0a1 - beta b0b1 = a0b0
	double unknown[2];
	double a0a1[3];
	double b0b1[3];
	double a0b0[3];

	for (int i = 0; i < 3; i++) {
		a0a1[i] = la1[i] - la0[i];
		b0b1[i] = lb1[i] - lb0[i];
		a0b0[i] = lb0[i] - la0[i];
	}

	double A[4];
	double a0a1_2D[2];
	double b0b1_2D[2];

	A[0] = a0a1_2D[0] = a0a1[(planeToProject + 1) % 3];
	A[2] = b0b1_2D[0] = -b0b1[(planeToProject + 1) % 3];
	unknown[0] = a0b0[(planeToProject + 1) % 3];
	A[1] = a0a1_2D[1] = a0a1[(planeToProject + 2) % 3];
	A[3] = b0b1_2D[1] = -b0b1[(planeToProject + 2) % 3];
	unknown[1] = a0b0[(planeToProject + 2) % 3];

	if (Parallel2D(a0a1_2D, b0b1_2D))
		return false;

	/*int dummy[2];
	lapack_int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 2, 1, A, 2, dummy, unknown, 2);
	if (info != 0) {
		cerr << "ERROR in MortarMath::intersect()! LAPACK error info: " << info << endl;
		exit(EXIT_FAILURE);
	}*/
	solve2x2LinearSystem(A, unknown);

	assert(unknown[1] > -EPS);
	assert(unknown[1] < 1+EPS);

	for (int i = 0; i < 3; i++)
		intersection[i] = la0[i] + unknown[0] * a0a1[i];
	return true;
}

bool MortarMath::Parallel2D(double *vec1, double *vec2) {
	const double EPS = 1E-10; // tolerance for parallel, if it is too small, intersection will be bullshit.
	double l1 = vectorLength2D(vec1);
	double l2 = vectorLength2D(vec2);
	for (int i=0; i<2; i++) {
		vec1[i] /= l1;
		vec2[i] /= l2;
	}
	double crossProduct = vec1[0] * vec2[1] - vec2[0] * vec1[1];
	return fabs(crossProduct) < EPS;
}

double MortarMath::vectorLength(const double *vec) {
	double length = 0;
	for (int i = 0; i < 3; i++)
		length += vec[i] * vec[i];
	length = sqrt(length);
	return length;
}

double MortarMath::vectorLength2D(const double *vec) {
	double length = 0;
	for (int i = 0; i < 2; i++)
		length += vec[i] * vec[i];
	length = sqrt(length);
	return length;
}

double MortarMath::dotProduct(const double *vec1, const double *vec2) {
	double product = 0;
	for (int i = 0; i < 3; i++)
		product += vec1[i] * vec2[i];
	return product;
}

int MortarMath::projectToPlane(const double *unitNormal) {
	double plane[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	double max = -1;
	int imax = -100;
	for (int i = 0; i < 3; i++) {
		double tmp = fabs(dotProduct(plane[i], unitNormal));
		if (tmp > max) {
			max = tmp;
			imax = i;
		}
	}
	return imax;
}

void MortarMath::solve2x2LinearSystem(const double *A, double *b) {
	// A is column major
	double detA = A[0]*A[3] - A[2]*A[1];
	double det0 = A[3]*b[0] - A[2]*b[1];
	double det1 = A[0]*b[1] - A[1]*b[0];
	b[0] = det0 / detA;
	b[1] = det1 / detA;
}

void MortarMath::solve3x3LinearSystem(const double *A, int planeToProject, double *b) {
	int l0 = (planeToProject + 1) % 3;
	int l1 = (planeToProject + 2) % 3;
	double detA1A2 = A[l0+3] * A[l1+6] - A[l1+3] * A[l0+6];
	double detA0A2 = A[l0] * A[l1+6] - A[l1] * A[l0+6];
	double detA0A1 = A[l0] * A[l1+3] - A[l1] * A[l0+3];
	double detA0b = A[l0] * b[l1] - A[l1] * b[l0];
	double detbA2 = b[l0] * A[l1+6] - b[l1] * A[l0+6];
	double detA1b = A[l0+3] * b[l1] - A[l1+3] * b[l0];
	double detA = detA1A2 - detA0A2 + detA0A1;
	double det0 = detA1A2 - detbA2 - detA1b;
	double det1 = detbA2 - detA0A2 + detA0b;
	double det2 = detA1b - detA0b + detA0A1;
	b[0] = det0 / detA;
	b[1] = det1 / detA;
	b[2] = det2 / detA;
}

void MortarMath::computeMatrixProduct(int n, const double *A, double *B) {
	double C[n*n];
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			double tmp = 0.0;
			for (int k=0; k<n; k++) {
				tmp += A[i*n + k] * B[j + k*n];
			}
			C[i*n + j] = tmp;
		}
	}
	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			B[i*n + j] = C[i*n + j];
}

void MortarMath::printElem(const double *elem, int size) {
	for (int i = 0; i < size; i++) {
		cout << endl;
		for (int j = 0; j < 3; j++)
			cout << elem[i * 3 + j] << "  ";
	}
	cout << endl;
}

void MortarMath::printPoint(const double *p) {
	cout << endl << "Coor: ";
	for (int i = 0; i < 3; i++) {
		cout << p[i] << "  ";
	}
	cout << endl;
}
