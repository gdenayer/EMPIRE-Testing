#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml > outputEmperor.log & EPID=$!
		sleep 1s

		blockMesh -case turekFSIBenchmarkOF
		pimpleDyMFsiFoam -case turekFSIBenchmarkOF >outputPimpleDyMFsiFoam.log &
		recorderClient > outputRecorderClient.log &

		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi