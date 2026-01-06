#!/bin/bash 

#SBATCH --job-name=lammps
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=27
#SBATCH --output=lammps_%j.out
#SBATCH --partition=lammps
#SBATCH --time=99999999:00:00

cd $SLURM_SUBMIT_DIR

#mpich
export PATH=$PATH:/opt/mpich3.3.2/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/mpich3.3.2/lib
export MANPATH=$MANPATH:/opt/mpich3.3.2/share/man
export PATH=$PATH:/opt/lammps-2Aug2023/bin

mpirun -np 27 lmp_intel_cpu_intelmpi -in in.input
