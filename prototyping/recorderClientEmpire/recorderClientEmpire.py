import os
from ctypes import cdll
from ctypes import *
import ctypes as c_types
import time
import numpy

# this tests recieve mesh calls of EMPIRE_API
print('------------------------------------')
print('|                                  |')
print('|    RecorderClientEmpire          |')
print('|                                  |')
print('|    Author : Andreas Apostolatos  |')
print('|                                  |')
print('------------------------------------\n')

# import EMPIRE_API
libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])

# Connect to EMPIRE using an xml file
libempire_api.EMPIRE_API_Connect("./recorderClientEmpire.xml")

# Read attributes of the client from the xml file
EMPIRE_API_getUserDefinedText = libempire_api.EMPIRE_API_getUserDefinedText;
EMPIRE_API_getUserDefinedText.restype = c_char_p
resName = EMPIRE_API_getUserDefinedText('resName');
resStep = EMPIRE_API_getUserDefinedText('resStep');
noTimeSteps = int(EMPIRE_API_getUserDefinedText('noTimeSteps'));
geomFile = EMPIRE_API_getUserDefinedText('geomFile');
resFile = EMPIRE_API_getUserDefinedText('resFile');

# String to find
stringToFind = 'Result "' + resName + '" "EMPIRE_CoSimulation" ' + resStep + ' vector OnNurbsSurface'
stringToStart = 'Values'
stingToEnd = 'End Values'

####
#### Reading an out.georhino.txt to get the number of Control Points per patch
####

## Define the out.georhino.txt file where the output geometry is defined
print 'Reading geometry file ', geomFile

## Initialize number of CPs in u- and v-directions
noUCPs = []
noVCPs = []

FoundLine = False
noTimes = 4
iTimes = 0

# Initialize counter
counter = 1

# Loop over all lines in the input file and store the number of CPs in u- and v-directions for each patch
with open(geomFile) as file:
	for line in file:
		if '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0' in line or FoundLine:
			FoundLine = True
			iTimes = iTimes + 1
		if FoundLine and iTimes == 4:
			# print line
			splitLine = line.split()
			noVCPs.append(int(splitLine[1]))
			noUCPs.append(int(splitLine[2]))
			iTimes = 0
			FoundLine = False
			counter += 1

# Number of patches
noPatches = len(noUCPs)
print '\tNumber of patches is ', noPatches

#####
##### Reading an out.post.res to get the values of the field
#####

# Define the out.post.res to get the values of the field
print 'Reading results file ', resFile

# Initialize field at each Cartesian direction
fieldX = []
fieldY = []
fieldZ = []

# Initialize number of patches
noPatches = 0;

# Loop over all lines of the file to get the field values at the chosen step
with open(resFile) as file:
	foundResults = False
	foundValues = False
	for line in file:
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
			if len(splitLine) > 1:
				fieldX.append(float(splitLine[0]))
				fieldY.append(float(splitLine[1]))
				fieldZ.append(float(splitLine[2]))
			else:
				noPatches = noPatches + 1
				continue
		else:
			AssertionError('No results found')
if foundResults:
	AssertionError('No results were found')
if foundValues:
	AssertionError('No values were found')

if len(fieldX) != len(fieldY) or len(fieldX) != len(fieldZ) or len(fieldY) != len(fieldZ):
	AssertionError('The results are not coherent')
noCPs = len(fieldX)
print '\tTotal number of interface CPs is ', noCPs

#####
##### Re-arranging the field values to comply with the numbering in Empire
#####
print 'Re-arranging the solution on the Control Points to match with the numbering in EMPIRE'

# Initialize the re-arranged vectors containing the field values 
fieldXRe = [None]*noCPs
fieldYRe = [None]*noCPs
fieldZRe = [None]*noCPs

# Initialize counter
counter = 0

# Initialize total number of CPs of previous patches
noCPsPrevious = 0

# Loop over all patches and subsequently over all Control Points of each patch to re-arrange the parametric directions u and v
for iPatches in range(noPatches):
	for iVPar in range(noVCPs[iPatches]):
		for iUPar in range(noUCPs[iPatches]):
			index = iUPar*noVCPs[iPatches] + iVPar
			fieldXRe[counter] = fieldX[noCPsPrevious + index]
			fieldYRe[counter] = fieldY[noCPsPrevious + index]
			fieldZRe[counter] = fieldZ[noCPsPrevious + index]
			counter = counter + 1
	noCPsPrevious += noUCPs[iPatches]*noVCPs[iPatches]

# Create the field to be sent
noDOFs = 3*noCPs
field = [None]*noDOFs
for i in range(len(fieldX)):
	field[3*i + 0] = float(fieldXRe[i])
	field[3*i + 1] = float(fieldYRe[i])
	field[3*i + 2] = float(fieldZRe[i])
	
# Create the field in a c_type to be send to EMPIRE
EMPIRE_Disp = (c_double*(noDOFs))(0.0)
for i in range(noDOFs):
	EMPIRE_Disp[i] = field[i]

####
#### Loop over all time steps and send field to EMPIRE
####
print 'Sending the field to EMPIRE'
for i in range(noTimeSteps):
	libempire_api.EMPIRE_API_sendDataField(" ", noDOFs,EMPIRE_Disp)

####
#### Disconcect
####
libempire_api.EMPIRE_API_Disconnect();