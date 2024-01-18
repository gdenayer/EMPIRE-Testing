/***********************************************************************************************//**
 * \file GiDFileIO.h
 * This file holds the class GiDFileIO
 * \author Tianyang Wang
 * \date 5/25/2012
 * \version alpha
 **************************************************************************************************/
#ifndef GIDFILEIO_H_
#define GIDFILEIO_H_

#include <string>

class GiDResultInventory;

/********//**
 * \brief Class GiDFileIO holds a number of static functions which are used to parse GiD format files
 ***********/
class GiDFileIO {
public:
    /********//**
     * \brief Read a .msh file, initialize mesh data
     *
     * \param[in] fileName name of the mesh file
     * \param[out] numNodes number of nodes
     * \param[out] numElems number of elements
     * \param[out] nodesPerElem number of nodes per element
     * \param[out] nodeCoors coordinates of nodes (x,y,z)
     * \param[out] nodeNumbers index/id of each node
     * \param[out] elemTable "element table" or "connectivity table"
     ***********/
    static void readDotMsh(std::string fileName, int &numNodes, int &numElems, int &nodesPerElem,
            double *&nodeCoors, int *&nodeNumbers, int *&elemTable);
    /********//**
     * \brief Write a .msh file
     *
     * \param[in] fileName name of the mesh file
     * \param[in] numNodes number of nodes
     * \param[in] numElems number of elements
     * \param[in] nodesPerElem number of nodes per element
     * \param[in] nodeCoors coordinates of nodes (x,y,z)
     * \param[in] nodeNumbers index/id of each node
     * \param[in] elemTable "element table" or "connectivity table"
     ***********/
    static void writeDotMsh(std::string fileName, int numNodes, int numElems, int nodesPerElem,
            double *nodeCoors, int *nodeNumbers, int *elemTable);

    /********//**
     * \brief Write a .res file, this function allows write several results to a .res file
     * in a single call, but the limitation is the analysis name and the step number is the
     * same among different results.
     *
     * This function should be deleted. Only testMapper in prototype calls it.(8/2/2012 by Tianyang)
     *
     * \param[in] fileName name of the .res file
     * \param[in] numNodes number of nodes
     * \param[in] nodeNumbers index/id of each node
     * \param[in] analysisName analysis name
     * \param[in] numResults number of results
     * \param[in] names array of result names
     * \param[in] results array of results
     ***********/
    static void writeDotRes(std::string fileName, int numNodes, int *nodeNumbers,
            std::string analysisName, int numResults, std::string *names, double **results);

    /********//**
     * \brief Initialize a .res instance. This function works together with the function
     * appendResultToDotRes(string,string,int,string,int,int*,double*) and writeDotRes(string)
     * \return a data structure representing a .res file
     ***********/
    static GiDResultInventory *initDotRes();

    /********//**
     * \brief Append result to the .res instance
     * \param[in] ri data structure representing a .res file
     * \param[in] resultName name of the result
     * \param[in] analysisName analysis name
     * \param[in] stepNum step number
     * \param[in] type "Vector" or "Scalar"
     * \param[in] numNodes number of nodes
     * \param[in] nodeNumbers index/id of each node
     * \param[in] data data of the result
     ***********/
    static void appendResultToDotRes(GiDResultInventory *ri, std::string resultName,
            std::string analysisName, int stepNum, std::string type, int numNodes, int *nodeNumbers,
            double *data);

    /********//**
     * \brief write the .res instance
     * \param[in] ri data structure representing a .res file
     * \param[in] fileName name of the .res file
     ***********/
    static void writeDotRes(GiDResultInventory *ri, std::string fileName);

    /********//**
     * \brief write error (can be not error) to a gnuplot file, this function has some limitations,
     * and has to be further improved. -> TODO
     *
     * \param[in] testCaseName name of the test case
     * \param[in] numNodes number of nodes
     * \param[in] nodeCoors coordinates of nodes (x,y,z)
     * \param[in] errorVec array of errors in a vector format (x,y,z)
     ***********/
    static void writeGnuplotError(std::string testCaseName, int numNodes, double *nodeCoors,
            double *errorVec);
};

#endif
