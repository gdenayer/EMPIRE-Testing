#include "FieldCreator.h"
#include <math.h>

using namespace std;

FieldCreator::FieldCreator(int _numNodes, int _numElems, int _nodesPerElem,
		double *_nodeCoors, int *_nodeNumbers, int *_elemTable,
		string _function) :
		AbstractDataCreator(_numNodes, _numElems, _nodesPerElem, _nodeCoors,
				_nodeNumbers, _elemTable,_function) {
	// empty
}

FieldCreator::~FieldCreator() {
	// empty
}

void FieldCreator::create(double *data) {
	//cout << " FieldCreator::create(double *data)" <<endl;
	for (int i = 0; i < numNodes; i++) {
	    data[i * DOF_DIM] = 0;
	    data[i * DOF_DIM + 1] = 0;
		double x = nodeCoors[i * COOR_DIM];
		double y = nodeCoors[i * COOR_DIM + 1];
		double z = nodeCoors[i * COOR_DIM + 2];
		data[i * DOF_DIM + 2] = funcZ(x, y, z);
	}
}
