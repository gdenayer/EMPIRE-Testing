/******************************************************************************//**
 * \file MortarMath.h
 * The header file of class MortarMath.
 * \author Tianyang Wang
 * \date 9/21/2011
 *********************************************************************************/
#ifndef MORTARMATH_H_
#define MORTARMATH_H_

class GaussQuadrature;
class Integrand;

/******************************************************************************//**
 * \brief math used by mortar including Gauss quadrature, local coordinates computation, intersection,
 * \brief small linear system solver and so on.
 *********************************************************************************/
class MortarMath {
public:

	/********//**
	 **********************************************************************************
	 * \brief Project a number of points to a plane
	 * \param[in] pointOnPlane a point on the plane
	 * \param[in] unitNormal unit normal of the plane
	 * \param[in] points the points to be projected
	 * \param[in] num number of points to be projected
	 * \param[in] projections the projections
	 ***********/
	static void projectToElemPlane(const double *pointOnPlane,
			const double *unitNormal, const double *points, int num,
			double *projections);
	/********//**
	 **********************************************************************************
	 * \brief Compute the normal vector of a triangle
	 * \param[in] triangle the triangle
	 * \param[in] unitLength whether set the length of the normal vector to 1 or not
	 * \param[out] normal normal vector
	 ***********/
	static void computeTriangleNormal(const double *triangle, bool unitLength,
			double *normal);
	/********//**
	 **********************************************************************************
	 * \brief Compute the area of a triangle
	 * \param[in] triangle the triangle
	 * \return the area
	 ***********/
	static double computeArea(const double *triangle);
	/********//**
	 **********************************************************************************
	 * \brief Compute local coordinates of a point in a triangle
	 * \param[in] triangle the triangle
	 * \param[in] planeToProject project to x-y or y-z or z-x plane (only use 2 coordinates)
	 * \param[in] point the point
	 * \param[out] localCoor local coordinates of the point
	 * \return a boolean saying whether the point is inside the triangle or not (true of false)
	 ***********/
	static bool computeLocalCoorInTriangle(const double *triangle, int planeToProject,
			const double *point, double *localCoor);
	/********//**
	 **********************************************************************************
	 * \brief Compute global coordinates of a point in a triangle
	 * \param[in] triangle the triangle
	 * \param[in] localCoor local coordinates of the point
	 * \param[out] globalCoor global coordinates of the point
	 ***********/
	static void computeGlobalCoorInTriangle(const double *triangle,
			const double *localCoor, double *globalCoor);
	/********//**
	 **********************************************************************************
	 * \brief Compute the shape function value by local coordinates in a quadrilateral
	 * \param[in] id the index/id of the shape function, could be 0, 1, 2 or 3
	 * \param[in] xi_eta local coordinates
	 * \return shape function value
	 ***********/
	static double getQuadShapeFunc(int id, const double *xi_eta);
	/********//**
	 **********************************************************************************
	 * \brief Compute the determinant of Jocobian matrix by local coordinates in a quadrilateral
	 * \param[in] quad the quadrilateral
	 * \param[in] xi_eta local coordinates
	 * \return determinant of Jocobian matrix
	 ***********/
	static double computeQuadDetJ(const double *quad, const double *xi_eta);
	/********//**
	 **********************************************************************************
	 * \brief Compute the normal vector of a quadrilateral
	 * \param[in] quad the quadrilateral
	 * \param[in] unitLength whether set the length of the normal vector to 1 or not
	 * \param[out] normal normal vector
	 ***********/
	static void computeQuadNormal(const double *quad, bool unitLength,
			double *normal);
	/********//**
	 **********************************************************************************
	 * \brief Compute the center of a quadrilateral
	 * \param[in] quad the quadrilateral
	 * \param[out] center the center
	 ***********/
	static void getQuadCenter(const double *quad, double *center);
	/********//**
	 **********************************************************************************
	 * \brief Compute local coordinates of a point in a quadrilateral
	 * \param[in] quad the quadrilateral
	 * \param[in] planeToProject project to x-y or y-z or z-x plane (only use 2 coordinates)
	 * \param[in] point the point
	 * \param[out] localCoor local coordinates of the point
	 * \return a boolean saying whether the point is inside the quadrilateral or not (true of false)
	 ***********/
	static bool computeLocalCoorInQuad(const double *quad, int planeToProject,
			const double *point, double *localCoor);

	/********//**
	 **********************************************************************************
	 * \brief Compute the square of the distance between two points
	 * \param[in] p1 the 1st point
	 * \param[in] p2 the 2nd point
	 * \return the square of the distance between p1 and p2
	 ***********/
	static double distanceSquare(const double *p1, const double *p2);
	/********//**
	 **********************************************************************************
	 * \brief Compute intersection between two lines
	 * \param[in] la1 the 1st point on la
	 * \param[in] la2 the 2nd point on la
	 * \param[in] lb1 the 1st point on lb
	 * \param[in] lb2 the 2nd point on lb
	 * \param[in] planeToProject project to x-y or y-z or z-x plane (only use 2 coordinates)
	 * \param[out] intersection the intersection of la and lb
	 * \return boolean saying whether the intersection is on lb or not (true of false)
	 ***********/
	static bool intersect(const double *la0, const double *la1,
			const double *lb0, const double *lb1, int planeToProject,
			double *intersection);
	/********//**
	 **********************************************************************************
	 * \brief Determine whether two vectors are parallel or not
	 * \param[in] vec1 the 1st vector (2D)
	 * \param[in] vec2 the 2nd vector (2D)
	 * \return boolean saying whether two vectors are parallel or not (true of false)
	 ***********/
	static bool Parallel2D(double *vec1, double *vec2);
	/********//**
	 **********************************************************************************
	 * \brief Compute the length of a vector
	 * \param[in] vec the vector
	 * \return the vector length
	 ***********/
	static double vectorLength(const double *vec);
	/********//**
	 **********************************************************************************
	 * \brief Compute the length of a 2D vector
	 * \param[in] vec the 2D vector
	 * \return the vector length
	 ***********/
	static double vectorLength2D(const double *vec);
	/********//**
	 **********************************************************************************
	 * \brief Compute the dot product of two vectors
	 * \param[in] vec1 the 1st vector
	 * \param[in] vec2 the 2nd vector
	 * \return the dot product
	 ***********/
	static double dotProduct(const double *vec1, const double *vec2);
	/********//**
	 **********************************************************************************
	 * \brief Project to x-y or y-z or z-x plane (only use 2 coordinates)
	 * \brief The result is the plane which has smallest angle with the unitNormal
	 * \param[in] unitNormal unit normal vector of a plane
	 * \return the id of the plane to be projected (0--->y-z, 1--->z-x, 2--->x-y)
	 ***********/
	static int projectToPlane(const double *unitNormal);
	/********//**
	 **********************************************************************************
	 * \brief Solve a 2x2 linear system by close form formula
	 * \param[in] A the left hand side matrix
	 * \param[in] b the right hand side vector
	 * \param[out] b the solution is written to b
	 ***********/
	static void solve2x2LinearSystem(const double *A, double *b);
	/********//**
	 **********************************************************************************
	 * \brief Solve a 3x3 linear system by close form formula, the system has one row which has all entries 1.
	 * \brief Therefore, it cannot solve general 3x3 linear system.
	 * \param[in] A the left hand side matrix
	 * \param[in] b the right hand side vector
	 * \param[in] planeToProject project to x-y or y-z or z-x plane (only use 2 coordinates)
	 * \param[out] b the solution is written to b
	 ***********/
	static void solve3x3LinearSystem(const double *A, int planeToProject, double *b);
	/********//**
	 **********************************************************************************
	 * \brief Compute the matrix product between two n*n matrices
	 * \param[in] n n*n matrices
	 * \param[in] A the 1st matrix
	 * \param[in] B the 2nd matrix
	 * \param[out] B B is overwritten by A*B
	 ***********/
	static void computeMatrixProduct(int n, const double *A, double *B);
	/********//**
	 **********************************************************************************
	 * \brief Print the coordinates of an element, used in debugging
	 * \param[in] elem the element
	 * \param[in] size 3---triangle, 4---quadrilateral
	 ***********/
	static void printElem(const double *elem, int size);
	/********//**
	 **********************************************************************************
	 * \brief Print the coordinates of a point, used in debugging
	 * \param[in] p the point
	 ***********/
	static void printPoint(const double *p);
};




/******************************************************************************//**
 * \brief Doing Gauss quadrature on a triangle or a quadrilateral
 *********************************************************************************/
class GaussQuadrature {
public:
    /********//**
     **********************************************************************************
     * \brief Do Gauss quadrature on a triangle
     * \param[in] integrand the integrand.
     * \return the integration value
     ***********/
    static double gaussQuadratureOnTriangle(Integrand *integrand);
    /********//**
     **********************************************************************************
     * \brief Do Gauss quadrature on a quadrilateral
     * \param[in] integrand the integrand.
     * \return the integration value
     ***********/
    static double gaussQuadratureOnQuad(Integrand *integrand);
    /********//**
     **********************************************************************************
     * \brief get the Gauss points (isoparametric coordinates) in a triangle
     * \param[in] numGaussPoints number of Gauss points
     * \return the Gauss points
     ***********/
    static const double* getTriangleGaussPoints(int numGaussPoints);
    /********//**
     **********************************************************************************
     * \brief get the weights of Gauss quadrature in a triangle
     * \param[in] numGaussPoints number of Gauss points
     * \return the weights
     ***********/
    static const double* getTriangleWeights(int numGaussPoints);
    /********//**
     **********************************************************************************
     * \brief get the weights of Gauss quadrature in a quadrilateral
     * \param[in] numGaussPoints number of Gauss points
     * \return the Gauss points
     ***********/
    static const double* getQuadGaussPoints(int numGaussPoints);
    /********//**
     **********************************************************************************
     * \brief get the weights of Gauss quadrature in a quadrilateral
     * \param[in] numGaussPoints number of Gauss points
     * \return the weights
     ***********/
    static const double* getQuadWeights(int numGaussPoints);
private:
    /// 3 Gauss points of a triangle
    static double* triGaussPoints3;
    /// weights on 3 Gauss points of a triangle
    static double* triWeights3;
    /// 6 Gauss points of a triangle
    static double* triGaussPoints6;
    /// weights on 6 Gauss points of a triangle
    static double* triWeights6;
    /// 7 Gauss points of a triangle
    static double* triGaussPoints7;
    /// weights on 7 Gauss points of a triangle
    static double* triWeights7;
    /// 12 Gauss points of a triangle
    static double* triGaussPoints12;
    /// weights on 12 Gauss points of a triangle
    static double* triWeights12;
    /// 1 Gauss point of a quadrilateral
    static double* quadGaussPoints1;
    /// weights on 1 Gauss point of a quadrilateral
    static double* quadWeights1;
    /// 4 Gauss points of a quadrilateral
    static double* quadGaussPoints4;
    /// weights on 4 Gauss points of a quadrilateral
    static double* quadWeights4;
    /// 9 Gauss points of a quadrilateral
    static double* quadGaussPoints9;
    /// weights on 9 Gauss points of a quadrilateral
    static double* quadWeights9;
};

/******************************************************************************//**
 * \brief The superclass (abstract) of all integrand.
 *********************************************************************************/
class Integrand {
private:
    /// number of Gauss points
    int numGaussPoints;
public:
    /********//**
     **********************************************************************************
     * \brief Destructor. There must be a virtual one so that the destructor of child classed
     * \brief is called! Otherwise, => memory leak!
     ***********/
    virtual ~Integrand() {
        // empty
    }
    /********//**
     **********************************************************************************
     * \brief Set to multiply which two shape functions (N? * N?),
     * \brief to be overloaded by integrand of multiplication of 2 shape functions.
     * \param[in] _shapeFunc1 index of the 1st shape function
     * \param[in] _shapeFunc2 index of the 2nd shape function
     ***********/
    virtual void setShapeFunctions(int shapeFunc1, int shapeFunc2) {
        // empty,
    }
    /********//**
     **********************************************************************************
     * \brief Evaluate the value of shape function product on a Gauss point
     * \param[in] gaussPointID id of the Gauss point.
     ***********/
    virtual double evaluate(int gaussPointID) = 0;
    /********//**
     **********************************************************************************
     * \brief Set the number of Gauss points
     * \param[in] num number of Gauss points
     ***********/
    void setNumGaussPoints(int num) {
        numGaussPoints = num;
    }
    /********//**
     **********************************************************************************
     * \brief Get the number of Gauss points
     * \return number of Gauss points
     ***********/
    int getNumGaussPoints() {
        return numGaussPoints;
    }
};

#endif /* MORTARMATH_H_ */
