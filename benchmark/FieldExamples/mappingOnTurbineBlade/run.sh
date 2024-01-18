#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		meshClientA meshClientAInput.xml >meshClientA.log &

		meshClientB meshClientBInput.xml >outputMeshClientB.log &
		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo   'Please use ICC.'
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi