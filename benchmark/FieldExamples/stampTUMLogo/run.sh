#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		dummyCSM dummyCSMInput >outputDummyCSM.log &
		blockMesh -case idealizedPistonTUM
		pimpleDyMFsiFoam -case idealizedPistonTUM >outputPimpleDyMFsiFoam.log &
		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		#mpiexec -np 1 ../../../EMPEROR/binICC/emperor emperorInput.xml & EPID=$!
		sleep 1s

		#mpiexec -np 1 ./../../clients/dummyCSM/binICC/dummyCSM dummyCSMInput >outputDummyCSM.log &

		#mpiexec -np 1 pimpleDyMFsiFoam  -case idealizedPistonTUM >outputPimpleDyMFsiFoam.log &
		
		#wait $EPID
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi