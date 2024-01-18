#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml > outputEmperor.log & EPID=$!
		sleep 1s

		blockMesh -case zeroThicknessMembraneOF
		pimpleDyMFsiFoam -case zeroThicknessMembraneOF > outputPimpleDyMFsiFoam.log &
		carat zeroThicknessMembraneCarat/zeroThicknessMembrane.dat >/dev/null 2>outputCarat.err &

		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi
