#ifndef GIDMESHINVENTORY_H_
#define GIDMESHINVENTORY_H_

#include <vector>
#include <ostream>
#include <string>

class GiDMesh;

using namespace std;

class GiDMeshInventory {
public:
	GiDMeshInventory(vector<GiDMesh*> *vecMesh);
	GiDMeshInventory(string meshFile);
	virtual ~GiDMeshInventory();
	vector<GiDMesh*> *getMeshes();

	void writeDotMsh(ostream &out);

private:
	vector<GiDMesh*> *_vecMesh;
	static string strToLower(string strInput);
};

#endif
