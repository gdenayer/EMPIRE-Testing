#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		#blockMesh -case OF
                pimpleDyMFsiOptFoam -case OF >/dev/null 2>outputPimpleDyMFsiOptFoam.err &
                carat carat/carat.dat >/dev/null 2>outputCarat.err &
                optimizationClient optClient.xml opt/out.georhino.txt 1>outputOpt.log 2>outputOpt.err &
		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi