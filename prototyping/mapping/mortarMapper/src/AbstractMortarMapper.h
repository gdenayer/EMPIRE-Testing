/******************************************************************************//**
 * \file AbstractMortarMapper.h
 * The header file of class AbstractMortarMapper.
 * \date 9/21/2011
 *********************************************************************************/
#ifndef ABSTRACTMORTARMAPPER_H_
#define ABSTRACTMORTARMAPPER_H_

#include <vector>
#include <set>
#include <map>

class ANNkd_tree;
class Integrand;

namespace flann {
    template<typename Distance> class Index;
    template<class T> struct L2;
    template<typename T> class Matrix;
}

/********//**
 * \brief The superclass of all types of mortar mappers.
 *********************************************************************************/
class AbstractMortarMapper {
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
	 * \param[in] _oppositeSurfaceNormal whether the interface of slave side and of master side have opposite normals  or not (true or false)
	 * \param[in] _dual whether or not to use dual mortar (true or false)
	 * \author Tianyang Wang
	 ***********/
	AbstractMortarMapper(int _slaveNumNodes, int _slaveNumElems,
			int _slaveNodesPerElem, double *_slaveNodeCoors,
			int *_slaveNodeNumbers, int *_slaveElemTable, int _masterNumNodes,
			int _masterNumElems, int _masterNodesPerElem,
			double *_masterNodeCoors, int *_masterNodeNumbers,
			int *_masterElemTable, bool _oppositeSurfaceNormal, bool _dual);
	/***********************************************************************************************
	 * \brief Destructor
	 * \author Tianyang Wang
	 ***********/
	virtual ~AbstractMortarMapper();
	/***********************************************************************************************
	 * \brief Do consistent mapping --- C_BB * masterDOF = C_BA * slaveDOF
	 * \param[in] slaveDOF the DOF of the slave side (e.g. x-displacements on all structure nodes)
	 * \param[out] masterDOF the DOF of the master side (e.g. x-displacements on all fluid nodes)
	 * \author Tianyang Wang
	 ***********/
	virtual void consistentMapping(double *slaveDOF, double *masterDOF);
    /***********************************************************************************************
     * \brief Do consistent mapping --- masterDOF = C_BA * slaveDOF
     * \param[in] slaveDOF the DOF of the slave side (e.g. x-tractions on all fluid nodes)
     * \param[out] masterDOF the DOF of the master side (e.g. x-forces on all structure nodes)
     * \author Tianyang Wang
     ***********/
    virtual void consistentMappingTraction2Force(double *slaveDOF, double *masterDOF);
    /***********************************************************************************************
	 * \brief Do conservative mapping --- slaveDOF = C_BA^T * C_BB^(-1) * masterDOF
	 * \param[in] masterDOF the DOF of the master side (e.g. x-forces on all fluid nodes)
	 * \param[out] slaveDOF the DOF of the slave side (e.g. x-forces on all structure nodes)
	 * \author Tianyang Wang
	 * ***********/
	virtual void conservativeMapping(double *masterDOF, double *slaveDOF);
    /***********************************************************************************************
     * \brief Do conservative mapping --- slaveDOF = C_BA^T * masterDOF
     * \param[in] masterDOF the DOF of the master side (REMARK: x-tractions on all fluid nodes)
     * \param[out] slaveDOF the DOF of the slave side (e.g. x-forces on all structure nodes)
     * \author Tianyang Wang
     * ***********/
    virtual void conservativeMappingTraction2Force(double *masterDOF, double *slaveDOF);
    /// defines number of threads used for MKL routines
    static int mklSetNumThreads;
    /// defines number of threads used for mapper routines
    static int mapperSetNumThreads;
protected:
	/***********************************************************************************************
	 * \brief This class is used to cooperate with the clipping algorithm
	 * ***********/
	class PolygonPointClipper {
	private:
		/// the polygon that clips other polygons
		const double *polygon;
		/// number of edges of the polygon
		int size;
		/// project to x-y or y-z or z-x plane (only use 2 coordinates)
		int planeToProject;
		/// flags used to determine inside/outside
		bool *insideFlag;
		/// slopes of all edges
		double *edgeSlopes;
		/// sign of whether to use y=f(x) or x=f(y)
		bool *reverseXY;

	public:
		/***********************************************************************************************
		 * \brief Constructor
		 * \param[in] _polygon the polygon that clips other polygons
		 * \param[in] _size number of nodes/edges in this polygon
		 * \param[in] _planeToProject project to x-y or y-z or z-x plane (only use 2 coordinates)
		 * \author Tianyang Wang
		 ***********/
		PolygonPointClipper(const double *_polygon, int _size,
				int _planeToProject);
		/***********************************************************************************************
		 * \brief Destructor
		 * \author Tianyang Wang
		 ***********/
		virtual ~PolygonPointClipper();
		/***********************************************************************************************
		 * \brief Decides whether a point is on the "inside" side of the i-th edge of the polygon
		 * \param[in] edgeID id of the edge
		 * \param[in] point the point
		 * \return if inside, return true, otherwise, return false
		 * \author Tianyang Wang
		 ***********/
		bool inside(int edgeID, double *point);
	};
	int slaveNumNodes;
	int slaveNumElems;
	int slaveNodesPerElem;
	double *slaveNodeCoors;
	int *slaveNodeNumbers;
	int *slaveElemTable;

	int masterNumNodes;
	int masterNumElems;
	int masterNodesPerElem;
	double *masterNodeCoors;
	int *masterNodeNumbers;
	int *masterElemTable;

	bool oppositeSurfaceNormal;
	bool dual;
	int N_P_E; // number of nodes per element

	/// directElemTable means the entries is not the node number, but the position in nodeCoors
	int *slaveDirectElemTable;
	int *masterDirectElemTable;
	/// given a node, all the elements containing it are got
	std::vector<int> **slaveNodeToElemTable;
	/// compute normals of all slave elements
	double *slaveElemNormals;


    flann::Index<flann::L2<double> > *FLANNkd_tree;
    flann::Matrix<double> *FLANNSlaveNodes;

	ANNkd_tree *slaveNodesTree;
	double **ANNSlaveNodes;

    /// C_BB csr format
    double *C_BB_A;
	/// C_BB csr format
	int *C_BB_IA;
	/// C_BB csr format
	int *C_BB_JA;
    /// dual version of C_BB_A, which is diagonal
    double *C_BB_A_DUAL;

	/// C_BA csr format
	double *C_BA_A;
	/// C_BA csr format
	int *C_BA_IA;
	/// C_BA csr format
	int *C_BA_JA;
	/// dual version of C_BA_A, it shares the same IA and JA with C_BA_A
	double *C_BA_A_DUAL;

	/// if teh meshes are triangulated, C_BB is factorC_BB=2.0 times larger
	double factorC_BB;
	/// if the meshes are triangulated, C_BA is factorC_BA=2.0*factorC_BB times larger
	double factorC_BA;

	/// pardiso variable
	void  *pt[64]; // this is related to internal memory management, see PARDISO manual
	/// pardiso variable
	int iparm[64];
	/// pardiso variable
	int mtype;
	/// pardiso variable
	int maxfct;
	/// pardiso variable
	int mnum;
	/// pardiso variable
	int msglvl;
	/// pardiso variable
	int neq;
	/// pardiso variable
	int nrhs;

	enum MeshLabel {
		SLAVE, MASTER
	};



	/***********************************************************************************************
	 * \brief Compute matrix C_BB
	 * \author Tianyang Wang
	 ***********/
	virtual void computeC_BB();
	/***********************************************************************************************
	 * \brief Compute matrix C_BA
	 * \author Tianyang Wang
	 ***********/
	virtual void computeC_BA();
	/***********************************************************************************************
	 * \brief Initialize pardiso to factorize C_BB
	 * \author Tianyang Wang
	 ***********/
	void initPardiso();
	/***********************************************************************************************
	 * \brief Deallocate the memory of pardiso
	 * \author Tianyang Wang
	 ***********/
	void deletePardiso();
	/***********************************************************************************************
	 * \brief Initialize all tables that help referring from an element to its nodes or vice versa
	 * \author Tianyang Wang
	 ***********/
	void initTables();
	/***********************************************************************************************
	 * \brief Deallocate the memory of all tables
	 * \author Tianyang Wang
	 ***********/
	void deleteTables();
	/***********************************************************************************************
	 * \brief Initialize the ANN nearest neighbor searching tree
	 * \author Tianyang Wang
	 ***********/
	void initANNTree();
	/***********************************************************************************************
	 * \brief Deallocate the memory of the searching tree
	 * \author Tianyang Wang
	 ***********/
	void deleteANNTree();
	/***********************************************************************************************
	 * \brief Set number of nodes per element. (assert both sides have the same value as n)
	 * \param[in] n set N_P_E = n
	 * \author Tianyang Wang
	 ***********/
	void setN_P_E(int n);

	/***********************************************************************************************
	 * \brief Clip the slaveElem by the masterElem. Do Gauss quadrature on the clipped.
	 * \param[in] masterElem the master element
	 * \param[in] planeToProject when do intersecting, project to x-y or y-z or z-x plane (only use 2 coordinates)
	 * \param[in] clipper the clipper constructed by the master element
	 * \param[in] slaveElem the slave element
	 * \param[out] result the Gauss quadrature of the shape function product. 9 entries if triangle, 16 entries if quadrilateral.
	 * \return if the two elements are overlapped, return true, otherwise, return false
	 * \author Tianyang Wang
	 ***********/
	bool clip(const double *masterElem, int planeToProject,
			AbstractMortarMapper::PolygonPointClipper &clipper,
			const double *slaveElem, double *result);
	/***********************************************************************************************
	 * \brief Compute the Gauss quadrature of the shape function product.
	 * \param[in] elem1 element 1
	 * \param[in] elem2 element 2
	 * \param[in] planeToProject when computing local coordinates, project to x-y or y-z or z-x plane (only use 2 coordinates)
	 * \param[in] polygon the clipped polygon
	 * \param[out] result the Gauss quadrature of the shape function product. 9 entries if triangle, 16 entries if quadrilateral.
	 * \author Tianyang Wang
	 ***********/
	void gaussQuadratureOnClip(const double *elem1, const double *elem2,
			int planeToProject, std::vector<double*> &polygon, double *result);
	/***********************************************************************************************
	 * \brief Compute the search radius to find the overlapping candidates of the master element.
	 * \param[in] masterElem
	 * \return the square of the search radius
	 * \author Tianyang Wang
	 ***********/
	double computeSearchRadiusSquare(const double* masterElem);
	/***********************************************************************************************
	 * \brief Find the overlapping candidates of the master element with the searching radius.
	 * \param[in] masterElem the master element
	 * \param[in] masterElemNormal the normal direction of the master element
	 * \param[in] searchRadiusSqr the square of the search radius
	 * \param[out] neighborElems the neighboring elements which are the overlapping candidates
	 * \author Tianyang Wang
	 ***********/
	void findCandidates(const double* masterElem,
			const double* masterElemNormal, double searchRadiusSqr,
			std::set<int> &neighborElems);
	/***********************************************************************************************
	 * \brief Kick out the candidates who have wrong normal direction. Only here the field "oppositeSurfaceNormal" is used.
	 * \brief This helps to avoid getting wrong overlap on the meshes like Turek benchmark 3D mesh.
	 * \param[in] masterUnitNormal the unit normal of the master element
	 * \param[in] slaveUnitNormal the unit normal of the slave element
	 * \param[in] bound the bound determining the upper limit of the normal direction difference
	 * \return should be kicked out or not (true or false).
	 * \author Tianyang Wang
	 ***********/
	bool kickOutCandidate(const double *masterUnitNormal,
			const double *slaveUnitNormal, double bound);
	/***********************************************************************************************
	 * \brief Project all nodes in the elements to the plane of the element.
	 * \param[in] elem the element (this can also be a point on the element plane)
	 * \param[in] planeUnitNormal the unit normal of the element plane
	 * \param[in] neighborElems the neighboring elements which are the overlapping candidates
	 * \param[out] projections the projections of all nodes in the elements
	 * \author Tianyang Wang
	 ***********/
	void projectToElemPlane(const double *elem, const double *planeUnitNormal,
	        std::set<int> &neighborElems, std::map<int, double*> &projections);

	/***********************************************************************************************
	 * \brief Given the element index/id, return the element
	 * \param[in] elemIndex the element index/id
	 * \param[in] label saying it is a slave element or a master element
	 * \param[out] elem the element corresponds to the element index
	 * \author Tianyang Wang
	 ***********/
	void getElemCoor(int elemIndex, MeshLabel label, double *elem);

	/***********************************************************************************************
	 * \brief Compute the table of normals of all slave elements. Must be implemented by concrete child classes.
	 * \author Tianyang Wang
	 ***********/
	virtual void computeSlaveElemNormals() = 0;
	/***********************************************************************************************
	 * \brief Compute shape function product on a single element. Must be implemented by concrete child classes.
	 * \param[in] elem the element
	 * \return the integrand for Gauss quadrature
	 * \author Tianyang Wang
	 ***********/
	virtual Integrand *getShapeFuncProductOnSingleElem(
			const double *elem) = 0;
	/***********************************************************************************************
	 * \brief Compute shape function product on two elements. Must be implemented by concrete child classes.
	 * \param[in] elem1 element 1
	 * \param[in] elem2 element 2
	 * \param[in] planeToProject when computing local coordinates, project to x-y or y-z or z-x plane (only use 2 coordinates)
	 * \param[in] clipTriangle a triangle of the clip (the clipped polygon is divided into triangles)
	 * \return the integrand for Gauss quadrature
	 * \author Tianyang Wang
	 ***********/
	virtual Integrand *getShapeFuncProductOn2Elems(
			const double *elem1, const double *elem2, int planeToProject,
			const double *clipTriangle) = 0;

	/***********************************************************************************************
	 * \brief Copy a point
	 * \param[in] origin to be copied
	 * \param[out] copy the copy
	 * \author Tianyang Wang
	 ***********/
	static void copyPoint(const double *origin, double *copy);
	/***********************************************************************************************
	 * \brief Copy a point
	 * \param[in] origin to be copied
	 * \param[out] copy the copy
	 * \author Tianyang Wang
	 ***********/
	static void copyElem(const double *origin, int size, double *copy);
	/***********************************************************************************************
	 * \brief Build a triangle
	 * \param[in] p0 the 1st point
	 * \param[in] p1 the 2nd point
	 * \param[in] p2 the 3rd point
	 * \param[out] triangle the triangle built by the three pointss
	 * \author Tianyang Wang
	 ***********/
	static void buildTrianagle(const double *p0, const double *p1,
			const double *p2, double *triangle);
	/***********************************************************************************************
	 * \brief Compute the center of a polygon
	 * \param[in] polygon the polygon
	 * \param[out] center the center of the polygon
	 * \author Tianyang Wang
	 ***********/
	static void polygonCenter(std::vector<double*> &polygon, double *center);
	/***********************************************************************************************
	 * \brief Compute the longest edge length of the element
	 * \param[in] elem the element
	 * \param[in] size 3--->triangle, 4--->quadrilateral
	 * \return the square of the longest length of the element
	 * \author Tianyang Wang
	 ***********/
	static double longestEdgeLengthSquare(const double *elem, int size);
	/***********************************************************************************************
	 * \brief Compute the coefficients to transfer the shape functions to dual shape functions
	 * \param[in] elem the element
	 * \param[out] coeffMatrix coefficient matrix
	 * \author Tianyang Wang
	 ***********/
	void computeDualCoeffMatrix(const double *elem, double *coeffMatrix);
    /***********************************************************************************************
     * \brief Check C_BB and C_BA related pointers, in order to verify the logic
     * \author Tianyang Wang
     ***********/
	void checkNullPointers();
};
#endif /* ABSTRACTMORTARMAPPER_H_ */
