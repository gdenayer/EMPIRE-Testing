#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml > outputEmperor.log & EPID=$!
		sleep 1s

		carat cbm_dynamic_nonlinear_newmark_shell8_FSI.txt >/dev/null 2>outputCarat.err &

		dummyCFD empireDummyCFD.xml >outputDummyCFD.log &
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		mpiexec -np 1 Emperor emperorInput.xml > outputEmperor.log & EPID=$!
		sleep 1s

		mpiexec -np 1 carat cbm_dynamic_nonlinear_newmark_shell8_FSI.txt >/dev/null 2>outputCarat.err &

		mpiexec -np 1 dummyCFD empireDummyCFD.xml >outputDummyCFD.log &
		
		wait $EPID
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi