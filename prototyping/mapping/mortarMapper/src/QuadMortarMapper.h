/***********************************************************************************************//**
 * \file QuadMortarMapper.h
 * The header file of class QuadMortarMapper.
 * \date 9/21/2011
 **************************************************************************************************/
#ifndef QUADMORTARMAPPER_H_
#define QUADMORTARMAPPER_H_

#include "AbstractMortarMapper.h"

/********//**
 * \brief Doing mortar mapping using quadrilateral shape functions (only worked on quadrilaterals).
 ***********/
class QuadMortarMapper: public AbstractMortarMapper {
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
	QuadMortarMapper(int _slaveNumNodes, int _slaveNumElems,
			int _slaveNodesPerElem, double *_slaveNodeCoors,
			int *_slaveNodeNumbers, int *_slaveElemTable, int _masterNumNodes,
			int _masterNumElems, int _masterNodesPerElem,
			double *_masterNodeCoors, int *_masterNodeNumbers,
			int *_masterElemTable, bool _oppositeSurfaceNormal, bool _dual);
	/***********************************************************************************************
	 * \brief Destructor
	 * \author Tianyang Wang
	 ***********/
	virtual ~QuadMortarMapper();

private:
	class ShapeFuncProductOnSameQuad;
	class ShapeFuncProductOn2Quads;

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

#endif /* QUADMORTARMAPPER_H_ */
