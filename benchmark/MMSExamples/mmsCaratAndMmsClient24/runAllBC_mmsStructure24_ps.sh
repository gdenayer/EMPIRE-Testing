#!/bin/bash


array_s=("0008" "0016" "0032" "0064" "0128" "0256")
#array_s=("0008")

#array_t=("0008" "0016" "0032" "0064" "0128")
array_t=("0001" "0002" "0004" "0008" "0016" "0032")
#array_t=("0008")

#array_tstep=("0.250000" "0.125000" "0.062500" "0.031250" "0.015625")
#array_tstep=("2.000" "1.000" "0.500" "0.250" "0.125")
#array_tstep=("0.250000")
rm -f probes_*

#done in the while loop now
#rm -f mmsClientInput.xml
#cp mmsClientInput.xml.orig mmsClientInput.xml

emperorFile='emperorInputBCprobes.xml'
caratBaseFile='carat_ps/rectFlat'
clientBaseFile='rectFlat'
mmsClientFile='mmsClientInput.xml'
pre='<nElemBC>'

echo "The array for spatial disretization contains ${#array_s[*]} entries"
echo "The array for temporal disretization contains ${#array_t[*]} entries"


j=0

while [ $j -lt ${#array_s[*]} ]  
do

rm -f mmsClientInput.xml emperorInputBCprobes.xml
cp mmsClientInput.xml.orig mmsClientInput.xml
cp emperorInputBCprobesIterative.xml.orig emperorInputBCprobes.xml
#cp emperorInputBCprobesLoose.xml.orig emperorInputBCprobes.xml

caratFile=$caratBaseFile${array_s[$j]}


###SPATIAL REFINEMENT
#target=$pre${array_s[$j]}
#target='<nElemBC>0008'
#sed -i s/'<nElemBC>0008'/$target/g $mmsClientFile

target=$clientBaseFile${array_s[$j]}
#target='rectFlat0128'
sed -i s/'rectFlat0008'/$target/g $mmsClientFile

sed -i s/'mmsStructureC'/'mmsStructure24'/g $mmsClientFile


###TIME REFINEMENT
timestep=`echo "scale=8; 1.0/${array_t[$j]}" | bc`
echo "The current time step size is ${timestep} entries"

a=${array_t[$j]} 
#b=1
#target=`echo "scale=0; $a+$b" | bc`
sed -i s/'10'/$a/g $emperorFile

a=${array_t[$j]} 
target='<numTimeSteps>'$a
#target=`echo "scale=0; $a+$b" | bc`
#target='<numTimeSteps>'$target
sed -i s/'<numTimeSteps>10'/$target/g $mmsClientFile

#target='<timeStep>'${array_tstep[$j]}
target='<timeStep>'$timestep
sed -i s/'<timeStep>0.00'/$target/g $mmsClientFile





	if [ $# -eq 1 ]; then
		if [ $1 = 'ICC' ]; then
			Emperor $emperorFile & EPID=$!
			sleep 1s

			carat $caratFile  >outputCarat.log &

			mmsClient $mmsClientFile >outputmmsClient.log &
		
			wait $EPID
	
		elif [ $1  = 'GCC' ]; then
			mpiexec -np 1 Emperor $emperorFile & EPID=$!
			sleep 1s

			carat $caratFile    >outputCarat.log &
		
			mpiexec -np 1 mmsClient $mmsClientFile >outputmmsClient.log &
		
			wait $EPID
	
		else
			echo   'Please input ICC or GCC.'
		fi
else
	echo   'Please input ICC or GCC.'
fi


mv probes probes_${array_s[$j]}





j=$[$j+1]
done



