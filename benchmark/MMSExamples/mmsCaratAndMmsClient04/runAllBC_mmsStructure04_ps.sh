#!/bin/bash


#array_s=("0008" "0016" "0032" "0064" "0128" "0256")
array_s=("0256")

rm -f probes_*

#done in the while loop now
#rm -f mmsClientInput.xml
#cp mmsClientInput.xml.orig mmsClientInput.xml

emperorFile='emperorInputBCprobesIterative.xml'
caratBaseFile='carat_ps/rectFlat'
clientBaseFile='rectFlat'
mmsClientFile='mmsClientInput.xml'
pre='<nElemBC>'

echo "The array for spatial disretization contains ${#array_s[*]} entries"



j=0

while [ $j -lt ${#array_s[*]} ]  
do

rm -f mmsClientInput.xml
cp mmsClientInput.xml.orig mmsClientInput.xml


#target=$pre${array_s[$j]}
#target='<nElemBC>0008'
#sed -i s/'<nElemBC>0008'/$target/g $mmsClientFile

target=$clientBaseFile${array_s[$j]}
#target='rectFlat0128'
sed -i s/'rectFlat0008'/$target/g $mmsClientFile

sed -i s/'mmsStructureC'/'mmsStructure04'/g $mmsClientFile


caratFile=$caratBaseFile${array_s[$j]}
#caratFile=carat/rectFlat0128
echo "carat input file has the source $caratFile"


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
#mv eProbes eProbes_${array_s[$j]}




j=$[$j+1]
done



