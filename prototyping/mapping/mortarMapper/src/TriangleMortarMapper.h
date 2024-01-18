/***********************************************************************************************//**
 * \file TriangleMortarMapper.h
 * The header file of class TriangleMortarMapper.
 * \date 9/21/2011
 **************************************************************************************************/
#ifndef TRIANGLEMORTARMAPPER_H_
#define TRIANGLEMORTARMAPPER_H_

#include "AbstractMortarMapper.h"

/********//**
 * \brief Doing mortar mapping using triangle shape functions (worked on both triangles and quadrilaterals).
 ***********/
class TriangleMortarMapper: public AbstractMortarMapper {
public:
	/***********************************************************************************************
	 * \brief Constructor
	 * \param[in] _slaveNumNodes number of nodes on slave side
	 * \param[in] _slaveNumElems number of elements on slave side
	 * \param[in] _slaveNodesPerElem number of nodes per element (only 3 or 4 is allowed)
	 * \param[in] _slaveNodeCoors x,y,z coordinates of all slave nodes
	 * \param[in] _slaveNodeNumbers the id of each slave node
	 * \param[in] _slaveElemTable the information of how the elements are constructed by the nodes
	 * \param[in] _master... analogical to _slave...
	 * \param[in] _oppositeSurfaceNormal whether the interface of slave side and master side have opposite normals
	 * \author Tianyang Wang
	 ***********/
	TriangleMortarMapper(int _slaveNumNodes, int _slaveNumElems, int _slaveNodesPerElem,
			double *_slaveNodeCoors, int *_slaveNodeNumbers,
			int *_slaveElemTable, int _masterNumNodes, int _masterNumElems,
			int _masterNodesPerElem, double *_masterNodeCoors,
			int *_masterNodeNumbers, int *_masterElemTable,
			bool _oppositeSurfaceNormal, bool _dual, bool symmetric);

	/***********************************************************************************************
	 * \brief Destructor
	 * \author Tianyang Wang
	 ***********/
	virtual ~TriangleMortarMapper();

private:
	class ShapeFuncProductOnSameTriangle;
	class ShapeFuncProductOn2Triangles;

	/***********************************************************************************************
	 * \brief Convert the quadrilateral mesh to triangular mesh.
	 * \param[in] if symmetric, divide a quadrilateral with both diagonals. If not, just take one.
	 * \author Tianyang Wang
	 ***********/
	void convertToTriangularMesh(bool symmetric); // convert the quadrilateral mesh to triangular mesh

	/***********************************************************************************************
	 * \brief Compute the table of normals of all slave elements. Implementing the abstract method of the super class.
	 * \author Tianyang Wang
	 ***********/
	void computeSlaveElemNormals();
	/***********************************************************************************************
	 * \brief Compute shape function product on a single element. Implementing the abstract method of the super class.
	 * \param[in] elem the element
	 * \return the integrand for Gauss quadrature
	 * \author Tianyang Wang
	 ***********/
	Integrand *getShapeFuncProductOnSingleElem(
				const double *elem);
	/***********************************************************************************************
	 * \brief Compute shape function product on two elements. Implementing the abstract method of the super class.
	 * \param[in] elem1 element 1
	 * \param[in] elem2 element 2
	 * \param[in] planeToProject when computing local coordinates, project to x-y or y-z or z-x plane (only use 2 coordinates)
	 * \param[in] clipTriangle a triangle of the clip (the clipped polygon is divided into triangles)
	 * \return the integrand for Gauss quadrature
	 * \author Tianyang Wang
	 ***********/
	Integrand *getShapeFuncProductOn2Elems(
				const double *elem1, const double *elem2, int planeToProject,
				const double *clipTriangle);
};

#endif
