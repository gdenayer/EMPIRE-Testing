#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml >outputEmperor.log & EPID=$!
		sleep 1s

		carat carat/cbm_beam1.dat 2>outputCarat.log &

		pimpleDyMFsiFoam -case OF >outputPimple.log &
		
		wait $EPID
	fi
else
	echo   'Please input ICC or GCC.'
fi
