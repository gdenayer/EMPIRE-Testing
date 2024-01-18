#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		blockMesh -case idealizedPiston
		decomposePar -force -case idealizedPiston
		mpirun -np 2 pimpleDyMFsiFoam -parallel -case idealizedPiston >outputPimpleDyMFsiFoam.err &
		dummyCSM dummyCSMInput >outputDummyCSM.log &

		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi