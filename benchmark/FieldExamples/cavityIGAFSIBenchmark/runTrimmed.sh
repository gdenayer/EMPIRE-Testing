#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml >outputEmperor.log & EPID=$!
		sleep 1s

		blockMesh -case cavityFSIBenchmarkOF
		pimpleDyMFsiFoam -case cavityFSIBenchmarkOF >/dev/null 2>outputPimpleDyMFsiFoam.err &
		carat cavityFSIBenchmarkCarat/cavityFSIBenchmarkCarat_trimmedSpline.dat >/dev/null 2>outputCarat.err &

		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi