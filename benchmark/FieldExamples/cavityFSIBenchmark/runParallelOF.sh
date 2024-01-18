#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		blockMesh -case cavityFSIBenchmarkOF
		decomposePar -force -case cavityFSIBenchmarkOF
		mpirun -np 3 pimpleDyMFsiFoam -parallel -case cavityFSIBenchmarkOF >/dev/null 2>outputPimpleDyMFsiFoam.err &
		carat cavityFSIBenchmarkCarat/cavityFSIBenchmarkCarat.dat >/dev/null 2>outputCarat.err &

		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi