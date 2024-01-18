from ctypes import *
import os
import math

# this tests recieve mesh calls of EMPIRE_API
print('------------------------------------')
print('|                                  |')
print('|    RecorderClientCarat           |')
print('|                                  |')
print('|    Author : Andreas Apostolatos  |')
print('|                                  |')
print('------------------------------------\n')

# import EMPIRE_API
libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
# libempire_api = cdll.LoadLibrary("absolute/path/to/libEMPIRE_API.so")# In case you don't want to use os module in python

EMPIRE_API_getUserDefinedText = libempire_api.EMPIRE_API_getUserDefinedText
EMPIRE_API_getUserDefinedText.restype = c_char_p
EMPIRE_API_sendMesh = libempire_api.EMPIRE_API_sendMesh
EMPIRE_API_sendMesh.argtypes = [POINTER(c_char_p), c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_int)]

# Connect to EMPIRE using an xml file
libempire_api.EMPIRE_API_Connect("./recorderClientCaratFEM.xml")

# Read attributes of the client from the xml file
resName = EMPIRE_API_getUserDefinedText('resName')
resStep = EMPIRE_API_getUserDefinedText('resStep')
noTimeSteps = int(EMPIRE_API_getUserDefinedText('noTimeSteps'))
noTimeStepsScaling = int(EMPIRE_API_getUserDefinedText('noTimeStepsScaling'))
geomFile = EMPIRE_API_getUserDefinedText('geomFile')
resFile = EMPIRE_API_getUserDefinedText('resFile')
nameCaratInput = EMPIRE_API_getUserDefinedText('nameCaratInput')
caratDirectory = EMPIRE_API_getUserDefinedText('caratDirectory')
phi = float(EMPIRE_API_getUserDefinedText('rotation'))

####
#### Read the mesh from the out.post.res file
####
print 'Reading mesh file ', geomFile

numNodes = 0
numElems = 0
nodeCoors = []
nodeIDs = []
numNodesPerElem = []
elemTable = []

# Initialize indicator flags
isCoordinates = False
isElements = False

with open(geomFile) as geometryFile:
	for line in geometryFile:
		if line.startswith('Coordinates'):
			isCoordinates = True
			continue
		if line.startswith('End Coordinates'):
			isCoordinates = False
		if isCoordinates:
			splitLine = line.split()
			numNodes += 1
			nodeIDs.append(int(splitLine[0]))
			for i in range(len(splitLine) - 1):
				nodeCoors.append(float(splitLine[i + 1]))
		if line.startswith('Elements'):
			isElements = True
			continue
		if line.startswith('End Elements'):
			isElements = False
		if isElements:
			splitLine = line.split()
			numElems += 1
			numNodesPerElem.append(len(splitLine) - 2)
			for i in range(len(splitLine) - 2):
				elemTable.append(int(splitLine[i + 1]))

####
#### Rotate the mesh nodes
####
if phi != 0:
	print 'Rotating mesh by ', phi, 'degrees'
	phiRad = float(phi)*float(math.pi)/float(180)
	for i in range(numNodes):
		XRotated = math.cos(phiRad)*nodeCoors[3*i + 0] + math.sin(phiRad)*nodeCoors[3*i + 1]
		YRotated = -math.sin(phiRad)*nodeCoors[3*i + 0] + math.cos(phiRad)*nodeCoors[3*i + 1]
		nodeCoors[3*i + 0] = XRotated
		nodeCoors[3*i + 1] = YRotated

####
#### Rotate the mesh nodes
####
print 'Sending mesh to Empire'

# Create the c type arrays
EMPIRE_numNodes = c_int(numNodes)
EMPIRE_numElems = c_int(numElems)
EMPIRE_nodeCoors = (c_double*(len(nodeCoors)))(0.0)
for i in range(len(nodeCoors)):
	EMPIRE_nodeCoors[i] = nodeCoors[i]
EMPIRE_nodeIDs = (c_int*(len(nodeIDs)))(0)
for i in range(len(nodeIDs)):
	EMPIRE_nodeIDs[i] = nodeIDs[i]
EMPIRE_numNodesPerElem = (c_int*(len(numNodesPerElem)))(0)
for i in range(len(numNodesPerElem)):
	EMPIRE_numNodesPerElem[i] = numNodesPerElem[i]
EMPIRE_elemTable = (c_int*(len(elemTable)))(0)
for i in range(len(elemTable)):
	EMPIRE_elemTable[i] = elemTable[i]

print EMPIRE_numNodes.value, len(EMPIRE_nodeCoors), len(EMPIRE_nodeIDs)
print EMPIRE_numElems.value, len(EMPIRE_numNodesPerElem), len(elemTable)
EMPIRE_meshName = c_char_p("myMesh1")
EMPIRE_API_sendMesh(EMPIRE_meshName, EMPIRE_numNodes, EMPIRE_numElems, EMPIRE_nodeCoors, EMPIRE_nodeIDs, EMPIRE_numNodesPerElem, EMPIRE_elemTable)

####
#### Reading an out.post.res to get the values of the field
####

# Define the out.post.res to get the values of the field
#resFile = './dataRecorderClient/out.post.res'
print 'Reading results file ', resFile

# String to find
stringToFind = 'Result "' + resName + '" "Eigenvector nmb." ' + resStep + ' Vector OnNodes'
stringToStart = 'Values'
stingToEnd = 'End Values'

# Initialize field at each Cartesian direction
fieldX = []
fieldY = []
fieldZ = []

# Loop over all lines of the file to get the field values at the chosen step
with open(resFile) as resultsFile:
	foundResults = False
	foundValues = False
	for line in resultsFile:
		if stringToFind in line:
			foundResults = True
			continue
		elif stingToEnd in line:
			foundResults = False
		
		if foundResults and not foundValues:
			if stringToStart in line:
				foundValues = True
				continue
			else:
				AssertionError('Results found but no string Values detected')
		elif foundResults and foundValues:
			splitLine = line.split()
			fieldX.append(float(splitLine[1]))
			fieldY.append(float(splitLine[2]))
			fieldZ.append(float(splitLine[3]))
		else:
			AssertionError('No results found')
if foundResults:
	AssertionError('No results were found')
if foundValues:
	AssertionError('No values were found')

if len(fieldX) != len(fieldY) or len(fieldX) != len(fieldZ) or len(fieldY) != len(fieldZ):
	AssertionError('The results are not coherent')
noCPs = len(fieldX)
print '\tTotal number of interface nodes is ', noCPs
			
# Create the field to be sent
noDOFs = 3*noCPs
field = [None]*noDOFs
normField = 0
for i in range(len(fieldX)):
	field[3*i + 0] = float(fieldX[i])
	field[3*i + 1] = float(fieldY[i])
	field[3*i + 2] = float(fieldZ[i])
	for j in range(3):
		normField += field[3*i + j]*field[3*i + j]

####
#### Rotating the mesh coordinates and the field
####
if phi != 0:
	print 'Rotating the field by ', phi, 'degrees'
	phiRad = float(phi)*float(math.pi)/float(180)
	for i in range(len(fieldX)):
		fieldX = math.cos(phiRad)*field[3*i + 0] + math.sin(phiRad)*field[3*i + 1]
		fieldY = -math.sin(phiRad)*field[3*i + 0] + math.cos(phiRad)*field[3*i + 1]
		field[3*i + 0] = fieldX
		field[3*i + 1] = fieldY

# Get the square root and invert the norm of the field
#normField = float(normField**(0.5))
#invNormField = float(1.0)/float(normField)

# Set a user specified initial scaling factor
#scalingInitial = float(1.0)/float(6.0)
#print 'Initial scaling factor equals ', scalingInitial

# Create the scaled field
#for i in range(noDOFs):
	#field[i] = invNormField*scalingInitial*field[i]

# Initialize the EMPIRE field
EMPIRE_Disp = (c_double*(noDOFs))(0.0)

####
#### Loop over all time steps and send field to EMPIRE
####
#print 'time step loop' 											#### Comment out for no time step looping ###
#for i in range(noTimeSteps): 									#### Comment out for no time step looping ###
	#print '\tTime step ', i										#### Comment out for no time step looping ###
	#scaling = float(i)/float(noTimeStepsScaling) 				#### Comment out for no time step looping ###
	#if scaling > 1: 											#### Comment out for no time step looping ###
		#scaling = 1 											#### Comment out for no time step looping ###
	#print '\tScaling field to be sent to EMPIRE by ', scaling	#### Comment out for no time step looping ###
scaling = 1 													#### Comment in for time step looping ###
for j in range(noDOFs): 										#### Comment out for no time step looping ###
	EMPIRE_Disp[j] = scaling* field[j]							#### Comment out for no time step looping ###
print '\tSending field to EMPIRE'
libempire_api.EMPIRE_API_sendDataField(" ", noDOFs, EMPIRE_Disp)

####
#### Disconcect
####
libempire_api.EMPIRE_API_Disconnect()
exit()