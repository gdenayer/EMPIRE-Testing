#ifndef GIDRESULTINVENTORY_H_
#define GIDRESULTINVENTORY_H_

#include <vector>
#include <ostream>
#include <strings.h>

class GiDResult;

using namespace std;

class GiDResultInventory {
public:
	GiDResultInventory(vector<GiDResult*> *vecResult);
	GiDResultInventory(string resFile);
	virtual ~GiDResultInventory();
	vector<GiDResult*> *getResults();
	void appendResults(GiDResultInventory *RI2);
	GiDResultInventory *mergeResults(GiDResultInventory *RI2);
	GiDResultInventory *copyResults();
	GiDResultInventory *expandScalarToVector();

	void writeDotRes(ostream &out);

private:
	vector<GiDResult*> *_vecResult;
	static string strToLower(string strInput);
};

#endif
