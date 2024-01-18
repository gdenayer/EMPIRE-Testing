#include <algorithm>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "GiDResultInventory.h"
#include "GiDBasicData.h"

GiDResultInventory::GiDResultInventory(vector<GiDResult*> *vecResult) :
		_vecResult(vecResult) {
}

GiDResultInventory::GiDResultInventory(string fileName) {
	ifstream resFile(fileName.c_str());
	if (!resFile) {
		cerr << ".res file not found" << '\n';
		exit(EXIT_FAILURE);
	}

	_vecResult = new vector<GiDResult*>();

	string dummy;
	string textLine;
	istringstream *lineStream = new istringstream("");
	while (getline(resFile, textLine, '\n')) {
		delete lineStream;
		lineStream = new istringstream(textLine.c_str());
		(*lineStream) >> dummy;
		/* parse the block "Result" **********************************************************************/

		if (strToLower(dummy).compare("result") == 0) {
			string name, analysisName, type, GPsName;
			int stepNum;
			bool onNodes;
			do {
				(*lineStream) >> dummy;
				name.append(" ").append(dummy);
			} while (name[name.size() - 1] != '\"'); // parse the name with spaces in
			name.erase(name.begin());
			do {
				(*lineStream) >> dummy;
				analysisName.append(" ").append(dummy);
			} while (analysisName[analysisName.size() - 1] != '\"'); // parse the analysisName with spaces in
			analysisName.erase(analysisName.begin());
			(*lineStream) >> stepNum;
			(*lineStream) >> type;
			(*lineStream) >> dummy;
			if (dummy.compare("OnNodes") == 0)
				onNodes = true;
			else
				onNodes = false;
			assert(onNodes);

			// start parsing the values
			do {
				getline(resFile, textLine, '\n');
				delete lineStream;
				lineStream = new istringstream(textLine.c_str());
				(*lineStream) >> dummy;
			} while (strToLower(dummy).compare("values") != 0); // find the line starting with keyword "values"
			int dim = GiDResult::getDOFDimByType(type);
			char e;
			vector<GiDValue*> *vecGiDValue = new vector<GiDValue*>();
			resFile >> e;
			while (tolower(e) != 'e') {
				resFile.putback(e);
				int id;
				resFile >> id;
				double* const value = new double[dim];
				for (int i = 0; i < dim; i++) {
					resFile >> value[i];
				}
				GiDValue * const mValue = new GiDValue(id, value);
				vecGiDValue->push_back(mValue);
				resFile >> e;
			}

			GiDResult * const result = new GiDResult(name, analysisName,
					stepNum, type, vecGiDValue);
			_vecResult->push_back(result);
		}
	}
}

GiDResultInventory::~GiDResultInventory() {
	for (unsigned i = 0; i < _vecResult->size(); i++)
		delete _vecResult->at(i);
	delete _vecResult;
}

vector<GiDResult*> *GiDResultInventory::getResults() {
	return _vecResult;
}

void GiDResultInventory::appendResults(GiDResultInventory *RI2) {
	GiDResult *result1 = this->getResults()->at(0);
	GiDResult *result2 = RI2->getResults()->at(0);
	assert(result1->getValues()->size() == result2->getValues()->size());
	vector<GiDResult*> *vecResult = this->getResults();
	for (unsigned i = 0; i < RI2->getResults()->size(); i++)
		vecResult->push_back(RI2->getResults()->at(i));
}

GiDResultInventory *GiDResultInventory::mergeResults(GiDResultInventory *RI2) {
	GiDResultInventory *RI10 = this->copyResults();
	GiDResultInventory *RI20 = RI2->copyResults();
	GiDResult *result1 = RI10->getResults()->at(0);
	GiDResult *result2 = RI20->getResults()->at(0);
	assert(result1->getValues()->size() == result2->getValues()->size());
	vector<GiDResult*> *vecResult = new vector<GiDResult*>;
	for (unsigned i = 0; i < RI10->getResults()->size(); i++)
		vecResult->push_back(RI10->getResults()->at(i));
	for (unsigned i = 0; i < RI20->getResults()->size(); i++)
		vecResult->push_back(RI20->getResults()->at(i));
	return new GiDResultInventory(vecResult);
}

GiDResultInventory *GiDResultInventory::copyResults() {
	vector<GiDResult*> *vecResult2 = new vector<GiDResult*>;
	vector<GiDResult*> *vecResult1 = this->getResults();

	for (unsigned i = 0; i < vecResult1->size(); i++) {
		string name = vecResult1->at(i)->getName();
		string analysisName = vecResult1->at(i)->getAnalysisName();
		int stepNum = vecResult1->at(i)->getStepNum();
		string type = vecResult1->at(i)->getType();
		vector<GiDValue*> *vecValue2 = new vector<GiDValue*>;
		vector<GiDValue*> *vecValue1 = vecResult1->at(i)->getValues();
		int dim = GiDResult::getDOFDimByType(type);
		for (unsigned j = 0; j < vecValue1->size(); j++) {
			int id = vecValue1->at(j)->getId();
			double *v = new double[dim];
			for (int k = 0; k < dim; k++)
				v[k] = (vecValue1->at(j)->getValue())[k];
			GiDValue *mvalue2 = new GiDValue(id, v);
			vecValue2->push_back(mvalue2);
		}
		GiDResult* result2 = new GiDResult(name, analysisName, stepNum, type,
				vecValue2);
		vecResult2->push_back(result2);
	}

	return new GiDResultInventory(vecResult2);
}

GiDResultInventory *GiDResultInventory::expandScalarToVector() {
	vector<GiDResult*> *vecResult2 = new vector<GiDResult*>;
	vector<GiDResult*> *vecResult1 = this->getResults();

	for (unsigned i = 0; i < vecResult1->size(); i++) {
		string name = vecResult1->at(i)->getName();
		name.insert(name.size() - 1, "(vector)");
		string analysisName = vecResult1->at(i)->getAnalysisName();
		int stepNum = vecResult1->at(i)->getStepNum();
		string type = vecResult1->at(i)->getType();
		vector<GiDValue*> *vecValue2 = new vector<GiDValue*>;
		vector<GiDValue*> *vecValue1 = vecResult1->at(i)->getValues();
		int dim = GiDResult::getDOFDimByType(type);
		if (dim != 1) {
			break;
		}

		type = "Vector";
		dim = GiDResult::getDOFDimByType(type);
		for (unsigned j = 0; j < vecValue1->size(); j++) {
			int id = vecValue1->at(j)->getId();
			double *v = new double[dim];

			v[0] = 0;
			v[1] = 0;
			v[2] = ((vecValue1->at(j)->getValue())[0] - 10.0) * 0.1;

			GiDValue *mvalue2 = new GiDValue(id, v);
			vecValue2->push_back(mvalue2);
		}
		GiDResult* result2 = new GiDResult(name, analysisName, stepNum, type,
				vecValue2);
		vecResult2->push_back(result2);
	}

	return new GiDResultInventory(vecResult2);
}

void GiDResultInventory::writeDotRes(ostream &out) {
	out << "GiD Post Results File 1.0\n\n";

	out << "#Edited by GiDParser:\n\n";

	vector<GiDResult*> *vecRes = getResults();
	for (unsigned i = 0; i < (*vecRes).size(); i++) {
		GiDResult *res = (*vecRes)[i];
		string dummy = "OnNodes";
		out << "Result\t" << res->getName() << '\t'
				<< res->getAnalysisName() << '\t' << res->getStepNum() << '\t'
				<< res->getType() << '\t' << dummy << '\n';
		out << "Values\n";
		int dim = GiDResult::getDOFDimByType(res->getType());
		vector<GiDValue*> *vecV = res->getValues();
		for (unsigned j = 0; j < (*vecV).size(); j++) {
			GiDValue *mv = (*vecV)[j];
			out << '\t' << mv->getId();
			for (int k = 0; k < dim; k++)
				out << '\t' << mv->getValue()[k];
			out << '\n';
		}
		out << "End Values\n\n";
	}
}

string GiDResultInventory::strToLower(string strInput) {
	string strLower;
	int i = 0;
	while (strInput[i]) {
		strLower.push_back(tolower(strInput[i]));
		i++;
	}
	return strLower;
}
