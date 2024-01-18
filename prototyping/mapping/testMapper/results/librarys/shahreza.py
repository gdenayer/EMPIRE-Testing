import os
import math
#	  Errors which can be used:
#	  raise AssertionError('No scalarfield in input-file')
#	  raise RuntimeError('Different meshsize, no mapping possible')
#	  raise ValueError('Different meshsize, no mapping possible')	



def run_command(command):
	#ATTENTION: every command will be started in a new shell -> FIX: combine more commands with ";"!
	os.system(command)
	print 'Running the command: ' + command
	
def run_script(script):
	os.system('sh ' + script)
	print 'Running this case is finished'
	
def publish_file(filePath):
	"""Publish the file to log.html"""
	print '*HTML* ' + '<a href=' + filePath + '> ' + filePath + ' </a>'

def plot_datafile(filePath):
	"""Call gnuplot to plot a data file"""
	"""Use gunplot to plot the file filepath.plot"""	
	os.system('gnuplot ' + filePath)

def publish_image(imagePath):
	"""Publish the image to log.html"""
	tmp = '<img src=' + imagePath + '>'
	print '*INFO* ' + tmp
	print '*HTML* ' + tmp
	print 'Image is published!'
	
	
def compare_pressure(path1, path2):
	pressure_field1 = scalarfield()
	pressure_field1.read_data_from_file(path1+'p')

	pressure_field2 = scalarfield()
	pressure_field2.read_data_from_file(path2+'p')

	#check if there is a different meshsize
	if pressure_field1.meshsize != pressure_field2.meshsize:
	  raise AssertionError('Different meshsize, no mapping possible')
	
	#create new scalarfield 'delta' for deviation between the pressure fields
	delta = scalarfield()
	for k in range(len(pressure_field1.data)):
	  delta.data.append(abs(pressure_field1.data[k] - pressure_field2.data[k]))
	print 'Maximum deviation in pressure: ' + str(delta.maximum_value())

	
def compare_velocity(path1, path2):
	velocity_field1 = vectorfield()
	velocity_field1.read_data_from_file(path1+'U')

	velocity_field2 = vectorfield()
	velocity_field2.read_data_from_file(path2+'U')	

	#check if there is a different meshsize
	if velocity_field1.meshsize != velocity_field2.meshsize:
	  raise AssertionError('Different meshsize, no mapping possible')

	#create new vectorfield 'delta' for deviation between the velocity fields
	delta = vectorfield()
	for k in range(len(velocity_field1.magnitude.data)):
	  delta.magnitude.data.append(abs(velocity_field1.magnitude.data[k] - velocity_field2.magnitude.data[k]))	  
	  
	print 'Maximum deviation in velocity magnitude: ' + str(delta.magnitude.maximum_value())
	print 'Mean deviation in velocity magnitude: ' + str(delta.magnitude.mean_value())