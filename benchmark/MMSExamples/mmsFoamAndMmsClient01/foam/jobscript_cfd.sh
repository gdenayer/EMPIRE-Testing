
#!/bin/sh

#PBS -N rectangular_g1_120k_0
#PBS -l nodes=2:ppn=2
###output files
#PBS -e rectangular_g1_120k_0_error
#PBS -o rectangular_g1_120k_0_out
###Mail to user
#PBS -m ae
#PBS -M fisch@bv.tum.de

cd $PBS_O_WORKDIR

mpiexec -n 4 movingPisoFoam -parallel &> output_rectangular_g1_120k_0

exit 0 