#!/bin/bash
#startFoam21


### VARIABLES & INPUT ###
emperorFile='emperorInputIterative.xml'
mmsClientFile='mmsClientInput.xml'
caratFile='carat/bottom'
clientBaseFile='bottom'
case=$(pwd)
#pre='<nElemBC>'

## Input Variables ##
mmsFunction='mmsFsi01'
totalTime=0.01

## List of spatial Discretizations (dominates the j-loop) ##
array_s=("0008" "0016" "0032" "0064" "0128" "0256" "0512")
#array_s=("0008")

## List of temporal Discretizations ##
#array_t=("0008" "0016" "0032" "0064" "0128")
array_t=("0001" "0002" "0004" "0008" "0016" "0032" "0064")
#array_t=("0001")

### INITIALIZING ###

cd ~/software/empire/build
make -j 3
cd ..
make -j 3 

cd ~/software/carat/src
make -j 3

cd /home/rupert/software/OpenFOAM/own-2.1.x/applications/solvers/mmsPimpleDyMFsiFoam
wmake


cd $case
rm -f probes* 
rm -rf  foam/probeset* foam/VTK*
echo "The array for spatial disretization contains ${#array_s[*]} entries"
echo "The array for temporal disretization contains ${#array_t[*]} entries"

if [ "${#array_s[*]}" != "${#array_t[*]}" ]
then
echo "DANGER! WRONG NUMBER OF ENTRIES IN THE ARRAYS FOR SPACE AND TIME"
exit -1
fi


### REFINEMENT LOOP STEPS ###
j=0
while [ $j -lt ${#array_s[*]} ]  
do

## Reset Input Files ##
rm -f mmsClientInput.xml emperorInputIterative.xml carat/bottom foam/system/controlDict 
rm -f foam/constant/polyMesh/blockMeshDict foam/constant/polyMesh/boundary foam/constant/polyMesh/faces foam/constant/polyMesh/neighbour foam/constant/polyMesh/owner foam/constant/polyMesh/points


cp mmsClientInputOrig.xml mmsClientInput.xml
cp emperorInputIterativeOrig.xml emperorInputIterative.xml
cp carat/bottomOrig carat/bottom
cp foam/constant/polyMesh/blockMeshDictOrig foam/constant/polyMesh/blockMeshDict
cp foam/system/controlDictOrig foam/system/controlDict

## Route carat Input File ##
#caratFile=$caratFile${array_s[$j]}
##^^ delete it

sed -i s/'mmsFunctionDummy'/$mmsFunction/g $mmsClientFile

## SPATIAL Refinement ##
s=${array_s[$j]}
target=$clientBaseFile$s
sed -i s/'meshFileDummy'/$target/g $mmsClientFile
sed -i s/'NELE'/$s/g $caratFile
sed -i s/'NELE'/$s/g foam/constant/polyMesh/blockMeshDict


## TEMPORAL Refinement ##
timestep=`echo "scale=8; $totalTime/${array_t[$j]}" | bc`
echo "The time step size is ${timestep}"

#a=${array_t[$j]} 
#b=1
#target=`echo "scale=0; $a+$b" | bc`
#sed -i s/'10'/$a/g $emperorFile

a=${array_t[$j]} 
sed -i s/'nTimeStepsDummy'/$a/g $mmsClientFile
sed -i s/'nTimeStepsDummy'/$a/g $emperorFile
sed -i s/'OutputTimestepDummy'/${array_t[$j]}/g foam/system/controlDict



sed -i s/'TimeStepDummy'/$timestep/g $mmsClientFile
sed -i s/'TimeStepDummy'/$timestep/g $caratFile
sed -i s/'TimeStepDummy'/$timestep/g foam/system/controlDict

sed -i s/'TotalTimeDummy'/$totalTime/g $caratFile
sed -i s/'TotalTimeDummy'/$totalTime/g foam/system/controlDict


    echo "--------------------------------------------------------------------------------"
    echo "Create Mesh for run case ${array_t[$j]}"
    echo "--------------------------------------------------------------------------------"
    rm -r foam/0
    cp -r foam/orig_0 foam/0
    blockMesh -case foam &> delete_me
    

	if [ $# -eq 1 ]; then
		if [ $1 = 'ICC' ]; then
		
    echo "--------------------------------------------------------------------------------"
    echo "Perform Calculation for run case ${array_t[$j]}"
    echo "--------------------------------------------------------------------------------"
		
			Emperor $emperorFile & EPID=$!
			sleep 1s

			carat $caratFile  >outputCarat.log &

			mmsFsiClient $mmsClientFile >outputmmsClient.log &

			mmsPimpleDyMFsiFoam -case foam >outputMmsPimpleDyMFsiFoam.log &
		    
		    
		wait $EPID
		
		
	echo "--------------------------------------------------------------------------------"
    echo "Build VTK files for run case ${array_t[$j]}"
    echo "--------------------------------------------------------------------------------"
		foamToVTK -case foam &> delete_me 
		
			#wait $EPID
	
		elif [ $1  = 'GCC' ]; then
			mpiexec -np 1 Emperor $emperorFile & EPID=$!
			sleep 1s

			carat $caratFile    >outputCarat.log &
		
			mpiexec -np 1 mmsFsiClient $mmsClientFile >outputmmsClient.log &
		
			blockMesh -case foam &
			mpiexec -np 1 mmsPimpleDyMFsiFoam -case foam >outputMmsPimpleDyMFsiFoam.log &
		    foamToVTK &> delete_me
			wait $EPID
	
		else
			echo   'Please input ICC or GCC.'
		fi
else
	echo   'Please input ICC or GCC.'
fi

    echo "--------------------------------------------------------------------------------"
    echo "Postprocessing for run case ${array_t[$j]}"
    echo "--------------------------------------------------------------------------------"
mv probesD probesD_${array_s[$j]}
mv probesE probesE_${array_s[$j]}
mv probesF probesF_${array_s[$j]}

mv $case/foam/probeset $case/foam/probeset_${array_s[$j]}
    
   			 cd $case/foam/probeset_${array_s[$j]}/0/
   			 if [ -f UDeviation ]
 			  then
   			   sed 's/(/ /g' UDeviation >UDev
    			  sed 's/)/ /g' UDev >UDeviation
    			  rm UDev
   			 fi
   			 cd $case/foam

   			 mv $case/foam/VTK $case/foam/VTK_${array_s[$j]}

  			  rm -r 1* 2* 3* 4* 5* 6* 7* 8* 9* &> delete_me
  			  rm delete_me

			 cd $case


j=$[$j+1]
done

#rm out*



