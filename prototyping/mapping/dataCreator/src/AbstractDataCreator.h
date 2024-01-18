/******************************************************************************//**
 * \file AbstractDataCreator.h
 * This file holds the class of the AbstractDataCreator
 * \author Tianyang Wang
 * \date 15/3/2012
 *********************************************************************************/

#ifndef ABSTRACTDATACREATOR_H_
#define ABSTRACTDATACREATOR_H_

/********//**
 * \brief This is the superclass of all solvers of the testAdapter
 ***********/

#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h>

class AbstractDataCreator {
protected:
    int numNodes;
    int numElems;
    int nodesPerElem;
    double *nodeCoors;
    int *nodeNumbers;
    int *elemTable;


    double (*funcX)(double, double, double);
    double (*funcY)(double, double, double);
    double (*funcZ)(double, double, double);
    static double sinFunction(double x, double y, double z) {
        return sin(6.285 * x) + sin(6.285 * y) + sin(6.285 * z);
    }
    static double constantFunction(double x, double y, double z) {
        return 1.0;
    }
    static double mms01X(double x, double y, double z) {
        return 1.0;
    }
    static double mms01Y(double x, double y, double z) {
        return 1.0;
    }
    static double mms01Z(double x, double y, double z) {
        return 1.0;
    }


public:
    /********//**
     **********************************************************************************
     * \brief Constructor
     ***********/
    AbstractDataCreator(int _numNodes, int _numElems, int _nodesPerElem, double *_nodeCoors,
            int *_nodeNumbers, int *_elemTable, std::string function) :
            numNodes(_numNodes), numElems(_numElems), nodesPerElem(_nodesPerElem), nodeCoors(
                    _nodeCoors), nodeNumbers(_nodeNumbers), elemTable(_elemTable) {


        if (function == "constant")
            funcZ = constantFunction;
        else if (function == "sine")
            funcZ = sinFunction;
        else if (function == "mms01"){
        	funcX = mms01X;
        	funcY = mms01Y;
        	funcZ = mms01Z;
        }

        else {
            std::cerr << std::endl
                    << "AbstractDataCreator::AbstractDataCreator: wrong function type!"
                    << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    /********//**
     **********************************************************************************
     * \brief Destructor
     ***********/
    virtual ~AbstractDataCreator() {
        //empty
    }

    /********//**
     **********************************************************************************
     * \brief Create data on a mesh
     * \param[out] data is the created (it is specified to be vector instead of scalar)
     ***********/
    virtual void create(double *data) = 0;
};

#endif /* ABSTRACTDATACREATOR_H_ */
