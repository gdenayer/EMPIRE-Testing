#include "GiDBasicData.h"
/* ------------------------- GiDNode -------------------------------*/
GiDNode::GiDNode(int id, double *coor) :
		_id(id), _coor(coor) {
}
GiDNode::~GiDNode() {
	delete[] _coor;
}
int GiDNode::getId() {
	return _id;
}
double* GiDNode::getCoordinate() {
	return _coor;
}

/* ------------------------- GiDElement -------------------------------*/
GiDElement::GiDElement(int id, int *nodeIds) :
		_id(id), _nodeIds(nodeIds) {
}
GiDElement::~GiDElement() {
	delete[] _nodeIds;
}
int GiDElement::getId() {
	return _id;
}
int *GiDElement::getNodeIds() {
	return _nodeIds;
}

/* ------------------------- GiDMesh -------------------------------*/
GiDMesh::GiDMesh(string name, string elementType, int dimension,
		int nodesPerElem, vector<GiDNode*> *nodes,
		vector<GiDElement*> *elements) :
		_name(name), _elementType(elementType), _dimension(dimension), _nodesPerElem(
				nodesPerElem), _nodes(nodes), _elements(elements) {
}
GiDMesh::~GiDMesh() {
	for (unsigned i = 0; i < _nodes->size(); i++)
		delete _nodes->at(i);
	for (unsigned i = 0; i < _elements->size(); i++)
		delete _elements->at(i);
	delete _nodes;
	delete _elements;
}
string GiDMesh::getName() {
	return _name;
}
string GiDMesh::getElementType() {
	return _elementType;
}
int GiDMesh::getDimension() {
	return _dimension;
}
int GiDMesh::getNodesPerElem() {
	return _nodesPerElem;
}
vector<GiDNode*> *GiDMesh::getNodes() {
	return _nodes;
}
vector<GiDElement*> *GiDMesh::getElements() {
	return _elements;
}
int GiDMesh::getNodeIndexByNodeId(int id) {
	if (id <= (int) _nodes->size())
		if (_nodes->at(id - 1)->getId() == id)
			return id - 1;
	for (int i = 0; i < (int) _nodes->size(); i++) {
		if (_nodes->at(i)->getId() == id)
			return i;
	}
	return -1;
}

/* ------------------------- GiDValue -------------------------------*/
GiDValue::GiDValue(int id, double *value) :
		_id(id), _value(value) {
}
GiDValue::~GiDValue() {
	delete[] _value;
}
int GiDValue::getId() {
	return _id;
}
double *GiDValue::getValue() {
	return _value;
}

/* ------------------------- GiDResult -------------------------------*/
GiDResult::GiDResult(string name, string analysisName, int stepNum, string type,
		vector<GiDValue*> *vecValue) :
		_name(name), _analysisName(analysisName), _stepNum(stepNum), _type(
				type), _vecValue(vecValue) {
	_DOFDim = getDOFDimByType(_type);
}
GiDResult::~GiDResult() {
	for (unsigned i = 0; i < _vecValue->size(); i++)
		delete _vecValue->at(i);
	delete _vecValue;
}

string GiDResult::getName() {
	return _name;
}
string GiDResult::getAnalysisName() {
	return _analysisName;
}
int GiDResult::getStepNum() {
	return _stepNum;
}
string GiDResult::getType() {
	return _type;
}

vector<GiDValue*> *GiDResult::getValues() {
	return _vecValue;
}

int GiDResult::getDOFDim() {
	return _DOFDim;
}

int GiDResult::getDOFDimByType(string type) {
	if (type == "Scalar")
		return 1;
	else if (type == "Vector")
		return 3;
	else {
		cerr << "GiDResult::GiDResult: Wrong result type!" << endl;
		exit(EXIT_FAILURE);
	}
}
