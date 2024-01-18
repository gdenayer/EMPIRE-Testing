#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		carat cavityFEMWithIGAInterfaceLayerCaratFEM/cavityFEMWithIGAInterfaceLayerCaratFEM_100Elem.dat >/dev/null 2>outputCaratFEM.err &
		blockMesh -case cavityFEMWithIGAInterfaceLayerOF
		pimpleDyMFsiFoam -case cavityFEMWithIGAInterfaceLayerOF >/dev/null 2>outputPimpleDyMFsiFoam.err &
# 		carat cavityFEMWithIGAInterfaceLayerCaratIGA/cavityFEMWithIGAInterfaceLayerCaratIGA.dat >/dev/null 2>outputCaratIGA.err &
		carat cavityFEMWithIGAInterfaceLayerCaratIGA/cavityFSIBenchmarkCarat_trimmedSpline.dat >/dev/null 2>outputCaratIGA.err &

		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi
