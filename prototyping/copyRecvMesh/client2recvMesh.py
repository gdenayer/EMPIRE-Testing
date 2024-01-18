import os
from ctypes import cdll
from ctypes import *
import ctypes as c_types
import time
import numpy


# this tests recieve mesh calls of EMPIRE_API
print "================================"
print "    This is copyRecv client!" 
print "================================"

# import EMPIRE_API
libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])

# import EMPIRE_API
libempire_api.EMPIRE_API_Connect("empireClient2recvMesh.xml")

# declaration of mesh initializers
# possible memory leakage due to memory allocation during EMPIRE_API_recvMesh
meshName=c_char_p("OPTMesh")
numNodes = c_int(0)
numElems = c_int(0)
nodes = POINTER(c_double)()
nodeIDs = POINTER(c_int)()
numNodesPerElem = POINTER(c_int)()
elems = POINTER(c_int)()


libempire_api.EMPIRE_API_recvMesh.argtypes = [POINTER(c_char_p), POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_double)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int))]
libempire_api.EMPIRE_API_recvMesh(meshName, byref(numNodes), byref(numElems), byref(nodes), byref(nodeIDs), byref(numNodesPerElem), byref(elems))

# print mesh
print "numNodes: ", numNodes
print "numElems: ", numElems
time.sleep(2)

print "nodes: "
for i in range(0,numNodes.value):
  print "Node: ", i+1, " x: ",nodes[i*3] ," y: ", nodes[i*3+1] ," z: ", nodes[i*3+2]
time.sleep(2)

print "nodeIDs: "  
for i in range(0,numNodes.value):
  print "Node: ", i+1, "  ", nodeIDs[i]
time.sleep(2)

count = 0
print "numNodesPerElem: "
for i in range(0,numElems.value):
  print "Elem: ", i+1, "  ", numNodesPerElem[i]
  count = count + numNodesPerElem[i]
time.sleep(2)

print "elems: "
for i in range(0,count):
  print "elems: ", i+1, "  ", elems[i]
  
libempire_api.EMPIRE_API_Disconnect()
exit()