#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		blockMesh -case OF
		pimpleDyMFsiOptFoam -case OF 1>outputOF.log 2>outputPimpleDyMFsiOptFoam.err &
		optimizationClient optClient.xml opt/out.georhino.txt 1>outputOpt.log 2>outputOpt.err &
		#carat cavityFSIBenchmarkCarat/cavityFSIBenchmarkCarat.dat >/dev/null 2>outputCarat.err &
		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi
