#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s
		
		pimpleDyMFsiFoam -case semisphericalIGAMembraneWithInternalPressureInFlowOF >/dev/null 2>outputPimpleDyMFsiFoam.err &
		carat semisphericalIGAMembraneWithInternalPressureInFlowCarat/semisphericalTrimmedIGAMembraneWithInternalPressureInFlow.dat >/dev/null 2>outputCarat.err &
		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi