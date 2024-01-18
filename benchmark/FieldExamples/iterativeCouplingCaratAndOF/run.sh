#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		carat idealizedPistonCarat/idealizedPistonCarat.dat >outputCarat.log &
		blockMesh -case idealizedPistonOF
		pimpleDyMFsiFoam -case idealizedPistonOF >outputPimpleDyMFsiFoam.log &
		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		#mpiexec -np 1 Emperor emperorInput.xml & EPID=$!
		sleep 1s

		#mpiexec -np 1 carat20.exe idealizedPistonCarat/idealizedPistonCarat.dat >outputCarat.log &
		#blockMesh -case idealizedPistonOF
		#mpiexec -np 1 pimpleDyMFsiFoam -case idealizedPistonOF >outputPimpleDyMFsiFoam.log &
		
		wait $EPID
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi