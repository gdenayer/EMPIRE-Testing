#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		dummyCSM dummyCSMInput >/dev/null 2>outputDummyCSM.err &
		blockMesh -case OF 
		pimpleDyMFsiFoam -case OF >/dev/null 2>outputPimpleDyMFsiFoam.err &
		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi