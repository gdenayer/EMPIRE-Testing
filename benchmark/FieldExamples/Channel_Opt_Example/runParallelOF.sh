#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s
		
		decomposePar -force -case channelOPTBenchmarkOF
		mpirun -np 2 adjointOptFoam -parallel -case channelOPTBenchmarkOF >OFoutput &
		carat channelOPTBenchmarkCarat/channelOPTBenchmarkCarat.dat >/dev/null 2>outputCarat.err &

		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi