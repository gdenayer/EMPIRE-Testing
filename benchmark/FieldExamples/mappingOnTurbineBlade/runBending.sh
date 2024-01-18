#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInputBending.xml & EPID=$!
		sleep 1s

		meshClientTurbomachinery meshClientTurbomachineryInput.xml > meshClientTurbomachinery.log &

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