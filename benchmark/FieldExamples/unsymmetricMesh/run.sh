#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		meshClientA meshClientAInput.xml >outputMeshClientA.log &

		meshClientB meshClientBInput.xml >outputMeshClientB.log &
		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		mpiexec -np 1 Emperor emperorInput.xml & EPID=$!
		sleep 1s

		mpiexec -np 1 meshClientA meshClientAInput.xml >outputMeshClientA.log &

		mpiexec -np 1 meshClientB meshClientBInput.xml >outputMeshClientB.log &
		
		wait $EPID
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi