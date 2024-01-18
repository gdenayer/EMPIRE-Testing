#!/bin/bash
if [ $# -eq 1 ]; then
	if [ $1 = 'ICC' ]; then
		#Emperor emperorInput.xml >outputEmperor.log & EPID=$!
		Emperor emperorInput.xml & EPID=$!
		sleep 1s

		carat semisphericalIGAInterfaceLayerCarat/semisphericalIGAInterfaceLayerCaratWithoutSingularities.dat >/dev/null 2>outputCaratIGA.err &

		pimpleDyMFsiFoam -case semisphericalFEMMembraneWithInternalPressureInFlowOF >/dev/null 2>outputPimpleDyMFsiFoam.err &
		
		carat semisphericalFEMMembraneWithInternalPressureInFlowCarat/semisphericalFEMMembraneWithInternalPressureInFlowCarat.dat >/dev/null 2>outputCaratFEM.err &
		#carat semisphericalFEMMembraneWithInternalPressureInFlowCarat/semisphericalFEMMembraneWithInternalPressureInFlowCaratFiner.dat >/dev/null 2>outputCaratFEM.err &
		#carat semisphericalFEMMembraneWithInternalPressureInFlowCarat/semisphericalFEMMembraneWithInternalPressureInFlowCaratFinest.dat >/dev/null 2>outputCaratFEM.err &
		
		wait $EPID
	
	elif [ $1  = 'GCC' ]; then
		echo "Stop openmp threading"
	
	else
		echo   'Please input ICC or GCC.'
	fi
else
	echo   'Please input ICC or GCC.'
fi
