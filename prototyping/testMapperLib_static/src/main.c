#include <stdio.h>
#include "MapperLib.h"

int main(){
    
  // mesh A data
  int AnumNodes = 6;
  int AnumElems = 2;
  int AnumNodesPerElem[2] = {4, 4};
  double Anodes[18] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    2.0, 0.0, 0.0,
    2.0, 1.0, 0.0
  };
  int AnodeIDs[6] = {1, 2, 3, 4, 5, 6};
  int Aelems[8] = {1, 2, 3, 4, 2, 5, 6, 3};
  
  // mesh B data
  int BnumNodes = 6;
  int BnumElems = 4;
  int BnumNodesPerElem[4] = {3, 3, 3, 3};
  double Bnodes[18] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    2.0, 0.0, 0.0,
    2.0, 1.0, 0.0
  };
  int BnodeIDs[6] = {1, 2, 3, 4, 5, 6};
  int Belems[12] = {1, 3, 4, 1, 2, 3, 2, 6, 3, 2, 5, 6};
  
  
  
  printf ("======================================================");
  printf ("\n----------------Initializing Mappers----------------"); 
  printf ("\n======================================================");
  printf ("\nInitializing NearestNeighborMapper");
  char NearestNeighborMapperName[] = "NearestNeighborMapper";
  init_FE_NearestNeighborMapper(NearestNeighborMapperName, AnumNodes, Anodes, BnumNodes, Bnodes);
  printf ("\nNearestNeighborMapper is initialized");
  
  printf ("\n----------------------------------------------------");
  printf ("\nInitializing NearestElementMapper");
  char NearestElementMapperName[] = "NearestElementMapper";
  init_FE_NearestElementMapper(NearestElementMapperName,
			       AnumNodes, AnumElems, AnumNodesPerElem, Anodes, AnodeIDs, Aelems, 
			       BnumNodes, BnumElems, BnumNodesPerElem, Bnodes, BnodeIDs, Belems);
  printf ("\nNearestElementMapper is initialized");
  
  printf ("\n----------------------------------------------------");
  printf ("\nInitializing BarycentricInterpolationMapper");
  char BarycentricInterpolationMapperName[] = "BarycentricInterpolationMapper";
  init_FE_BarycentricInterpolationMapper(BarycentricInterpolationMapperName, 
					 AnumNodes, Anodes, BnumNodes, Bnodes);
  printf ("\nBarycentricInterpolationMapper is initialized");
  
  printf ("\n----------------------------------------------------");
  printf ("\nInitializing MortarMapper");
  // mortar mapper options
  int oppositeSurfaceNormal = 0;
  int dual = 0;
  int enforceConsistency = 0;
  char MortarMapperName[] = "MortarMapper";
  init_FE_MortarMapper(MortarMapperName,
		       AnumNodes, AnumElems, AnumNodesPerElem, Anodes, AnodeIDs, Aelems, BnumNodes, 
		       BnumElems, BnumNodesPerElem, Bnodes, BnodeIDs, Belems, oppositeSurfaceNormal, 
		       dual, enforceConsistency);
  printf ("\nMortarMapper is initialized");
  
  printf ("\n======================================================");
  printf ("\n-----------------Running Mapper Tests-----------------"); 
  printf ("\n======================================================");
  printf ("\nInitializing data fields") ;
  int dimension = 1;
  const int dataSizeA = AnumNodes;
  const int dataSizeB = BnumNodes;
  const int dataSizeC = dataSizeA;
  double dataA[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  printf("\ndataA:");
  for (int i=0; i<dataSizeA; i++)    printf(" %E",dataA[i]);
  printf ("\n------------------------------------------------------");
  printf ("\nPerforming NearestNeighbor Mapping");
  double dataB_NN[6];
  double dataC_NN[6];
  
  // consistent mapping dataA -> dataB_NN
  printf ("\ndoConsistentMapping dataA->dataB_NN");
  doConsistentMapping(NearestNeighborMapperName, dimension, dataSizeA, dataA, dataSizeB, dataB_NN);
  printf("\ndataB_NN:");
  for (int i=0; i<dataSizeB; i++)    printf(" %E",dataB_NN[i]);
  
  // conservative mapping dataB_NN -> dataC_NN
  printf ("\ndoConservativeMapping dataB_NN->dataC_NN");
  doConservativeMapping(NearestNeighborMapperName, dimension, dataSizeB, dataB_NN, dataSizeC, dataC_NN);
  printf("\ndataC_NN:");
  for (int i=0; i<dataSizeC; i++)    printf(" %E",dataC_NN[i]);
  
  printf ("\n------------------------------------------------------");
  printf ("\nPerforming NearestElement Mapping");
  double dataB_NE[6];
  double dataC_NE[6];
  
  // consistent mapping dataA -> dataB_NE
  printf ("\ndoConsistentMapping dataA->dataB_NE");
  doConsistentMapping(NearestElementMapperName, dimension, dataSizeA, dataA, dataSizeB, dataB_NE);
  printf("\ndataB_NE:");
  for (int i=0; i<dataSizeB; i++)    printf(" %E",dataB_NE[i]);
  
  // conservative mapping dataB_NE -> dataC_NE
  printf ("\ndoConservativeMapping dataB_NE->dataC_NE");
  doConservativeMapping(NearestElementMapperName, dimension, dataSizeB, dataB_NE, dataSizeC, dataC_NE);
  printf("\ndataC_NE:");
  for (int i=0; i<dataSizeC; i++)    printf(" %E",dataC_NE[i]);
  
  printf ("\n------------------------------------------------------");
  printf ("\nPerforming BarycentricInterpolation Mapping");
  double dataB_BCI[6];
  double dataC_BCI[6];
  
  // consistent mapping dataA -> dataB_BCI
  printf ("\ndoConsistentMapping dataA->dataB_BCI");
  doConsistentMapping(NearestElementMapperName, dimension, dataSizeA, dataA, dataSizeB, dataB_BCI);
  printf("\ndataB_BCI:");
  for (int i=0; i<dataSizeB; i++)    printf(" %E",dataB_BCI[i]);
  
  // conservative mapping dataB_BCI -> dataC_BCI
  printf ("\ndoConservativeMapping dataB_BCI->dataC_BCI");
  doConservativeMapping(NearestElementMapperName, dimension, dataSizeB, dataB_BCI, dataSizeC, dataC_BCI);
  printf("\ndataC_BCI:");
  for (int i=0; i<dataSizeC; i++)    printf(" %E",dataC_BCI[i]);
 
  printf ("\n------------------------------------------------------");
  printf ("\nPerforming Mortar Mapping");
  double dataB_Mortar[6];
  double dataC_Mortar[6];
  
  // consistent mapping dataA -> dataB_Mortar
  printf ("\ndoConsistentMapping dataA->dataB_Mortar");
  doConsistentMapping(MortarMapperName, dimension, dataSizeA, dataA, dataSizeB, dataB_Mortar);
  printf("\ndataB_Mortar:");
  for (int i=0; i<dataSizeB; i++)    printf(" %E",dataB_Mortar[i]);
  
  // conservative mapping dataB_mortar -> dataC_mortar
  printf ("\ndoConservativeMapping dataB_Mortar->dataC_Mortar");
  doConservativeMapping(MortarMapperName, dimension, dataSizeB, dataB_Mortar, dataSizeC, dataC_Mortar);
  printf("\ndataC_Mortar:");
  for (int i=0; i<dataSizeC; i++)    printf(" %E",dataC_Mortar[i]);

  printf ("\n======================================================\n");
  
  // delete all the mappers
  deleteAllMappers();
  
  return 0;
  
}
