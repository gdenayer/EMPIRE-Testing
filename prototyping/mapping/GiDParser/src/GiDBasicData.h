#ifndef GIDBASICDATA_H_
#define GIDBASICDATA_H_

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>

using namespace std;

class GiDNode {
public:
	GiDNode(int id, double *coor);
	virtual ~GiDNode();
	int getId();
	double* getCoordinate();

private:
	int _id;
	double *_coor;
};

class GiDElement {
public:
	GiDElement(int id, int *nodeIds);
	virtual ~GiDElement();
	int getId();
	int *getNodeIds();

private:
	int _id;
	int *_nodeIds;
};

class GiDMesh {
public:
	GiDMesh(string name, string elementType, int dimension, int nodesPerElem,
			vector<GiDNode*> *nodes, vector<GiDElement*> *elements);
	virtual ~GiDMesh();
	string getName();
	string getElementType();
	int getDimension();
	int getNodesPerElem();
	vector<GiDNode*> *getNodes();
	vector<GiDElement*> *getElements();
	int getNodeIndexByNodeId(int id);

private:
	string _name;
	string _elementType;
	int _dimension;
	int _nodesPerElem;
	vector<GiDNode*> *_nodes;
	vector<GiDElement*> *_elements;
};

class GiDValue {
public:
	GiDValue(int id, double *value);
	virtual ~GiDValue();
	int getId();
	double *getValue();

private:
	int _id; // id of node
	double *_value;
};

class GiDResult {
public:
	GiDResult(string name, string analysisName, int stepNum, string type,
			vector<GiDValue*> *vecValue);
	virtual ~GiDResult();
	string getName();
	string getAnalysisName();
	int getStepNum();
	string getType();
	int getDOFDim();
	static int getDOFDimByType(string type);
	vector<GiDValue*> *getValues();

private:
	string _name;
	string _analysisName;
	int _stepNum;
	string _type;
	int _DOFDim;
	vector<GiDValue*> *_vecValue;
};

#endif
