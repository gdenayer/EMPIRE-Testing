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

# Connect to EMPIRE using an xml file
libempire_api.EMPIRE_API_Connect("./recorderClientCarat.xml")

# Read attributes of the client from the xml file
EMPIRE_API_getUserDefinedText = libempire_api.EMPIRE_API_getUserDefinedText
EMPIRE_API_getUserDefinedText.restype = c_char_p
resName = EMPIRE_API_getUserDefinedText('resName')
resStep = EMPIRE_API_getUserDefinedText('resStep')
noTimeSteps = int(EMPIRE_API_getUserDefinedText('noTimeSteps'))
noTimeStepsScaling = int(EMPIRE_API_getUserDefinedText('noTimeStepsScaling'))
geomFile = EMPIRE_API_getUserDefinedText('geomFile')
resFile = EMPIRE_API_getUserDefinedText('resFile')
nameCaratInput = EMPIRE_API_getUserDefinedText('nameCaratInput')
caratDirectory = EMPIRE_API_getUserDefinedText('caratDirectory')
phi = float(EMPIRE_API_getUserDefinedText('rotation'))

# String to find
stringToFind = 'Result "' + resName + '" "Eigenvector nmb." ' + resStep + ' Vector OnNodes'
stringToStart = 'Values'
stingToEnd = 'End Values'

####
#### Reading an out.georhino.txt to get the number of Control Points per patch
####

# Define the out.georhino.txt file where the output geometry is defined
# geomFile = './dataRecorderClient/out.georhino.txt'
print 'Reading geometry file ', geomFile

# Initialize number of CPs in u- and v-directions
noUCPs = []
noVCPs = []

# Loop over all lines in the input file and store the number of CPs in u- and v-directions for each patch
with open(geomFile) as file:
	for line in file:
		if 'MCTRL' in line:
			splitLine = line.split()
			splitPrevLine = prevLine.split()
			noVCPs.append(int(splitLine[2]) + 1)
			noUCPs.append(int(splitPrevLine[2]) + 1)
		prevLine = line

# Number of patches
noPatches = len(noUCPs)

####
#### Reading an out.post.res to get the values of the field
####

# Define the out.post.res to get the values of the field
#resFile = './dataRecorderClient/out.post.res'
print 'Reading results file ', resFile

# Initialize field at each Cartesian direction
fieldX = []
fieldY = []
fieldZ = []

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
print '\tTotal number of interface CPs is ', noCPs

####
#### Re-arranging the field values to comply with the numbering in Empire
####
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
	for iUPar in range(noUCPs[iPatches]):
		for iVPar in range(noVCPs[iPatches]):
			# print(iVPar*noUCPs[iPatches] + iUPar)
			# print(counter)
			index = iVPar*noVCPs[iPatches] + iUPar
			fieldXRe[counter] = fieldX[noCPsPrevious + index]
			fieldYRe[counter] = fieldY[noCPsPrevious + index]
			fieldZRe[counter] = fieldZ[noCPsPrevious + index]
			counter = counter + 1
	noCPsPrevious = noUCPs[iPatches]*noVCPs[iPatches]
			
# Create the field to be sent
noDOFs = 3*noCPs
field = [None]*noDOFs
normField = 0
for i in range(len(fieldX)):
	field[3*i + 0] = float(fieldXRe[i])
	field[3*i + 1] = float(fieldYRe[i])
	field[3*i + 2] = float(fieldZRe[i])
	for j in range(3):
		normField += field[3*i + j]*field[3*i + j]

####
#### Rotating the field
####
if phi != 0:
	print 'Rotating the field by ', phi, 'degrees'
	phi = float(phi)*float(math.pi)/float(180)
	for i in range(len(fieldX)):
		fieldX = math.cos(phi)*field[3*i + 0] + math.sin(phi)*field[3*i + 1]
		fieldY = -math.sin(phi)*field[3*i + 0] + math.cos(phi)*field[3*i + 1]
		field[3*i + 0] = fieldX
		field[3*i + 1] = fieldY
		
	print 'Rotating the geometry by ', phi, 'degrees'
	inputFilePrefix = caratDirectory+nameCaratInput
	inputFileName = inputFilePrefix + '.txt'
	isExistent = os.path.isfile(inputFileName)
	if not isExistent:
		AssertionError('File ', inputFileName, ' could not be found')
	inputFileRotated = inputFilePrefix + 'Rotated' + '.txt'
	
	inputFile = open(inputFileName, "r")
	outputFile = open(inputFileRotated, "w")
	
	output_lines = []
	for line in inputFile:
		if line.startswith(' CTRLPT  '):
			# Split line
			splitLine = line.split()
			
			# read the coordinates of the Control Point
			xCoord = float(splitLine[2])
			yCoord = float(splitLine[3])
			zCoord = float(splitLine[4])
			
			# Rotate them
			xCoordRot = math.cos(phi)*xCoord + math.sin(phi)*yCoord
			yCoordRot = -math.sin(phi)*xCoord + math.cos(phi)*yCoord
			zcoordRot = zCoord
			
			# Write the rotated values in the line
			lineRotated = ' '+splitLine[0]+'  '+splitLine[1]+'  '+str(xCoordRot)+'  '+str(yCoordRot)+'  '+str(zcoordRot)+'  '+splitLine[5] + '\n'
			output_lines.append(lineRotated)
			
		else:
			output_lines.append(line)
	print 'Input file with rotated Control Points written in ', inputFileRotated
	outputFile.writelines(output_lines)
	
	# Close the files
	inputFile.close()
	outputFile.close()

# Get the square root and invert the norm of the field
normField = float(normField**(0.5))
invNormField = float(1.0)/float(normField)

# Set a user specified initial scaling factor
scalingInitial = float(1.0)/float(6.0)
print 'Initial scaling factor equals ', scalingInitial

# Create the scaled field
for i in range(noDOFs):
	field[i] = invNormField*scalingInitial*field[i]

# Initialize the EMPIRE field
EMPIRE_Disp = (c_double*(noDOFs))(0.0)

####
#### Loop over all time steps and send field to EMPIRE
####
print 'time step loop' 										#### Comment out for no time step looping ###
for i in range(noTimeSteps): 									#### Comment out for no time step looping ###
	print '\tTime step ', i									#### Comment out for no time step looping ###
	scaling = float(i)/float(noTimeStepsScaling) 				#### Comment out for no time step looping ###
	if scaling > 1: 											#### Comment out for no time step looping ###
		scaling = 1 											#### Comment out for no time step looping ###
	print '\tScaling field to be sent to EMPIRE by ', scaling	#### Comment out for no time step looping ###
# scaling = 1 													#### Comment in for time step looping ###
for j in range(noDOFs): 									#### Comment out for no time step looping ###
	EMPIRE_Disp[j] = scaling* field[j]						#### Comment out for no time step looping ###
print '\tSending field to EMPIRE'
libempire_api.EMPIRE_API_sendDataField(" ", noDOFs, EMPIRE_Disp)

####
#### Disconcect
####
libempire_api.EMPIRE_API_Disconnect()
exit()