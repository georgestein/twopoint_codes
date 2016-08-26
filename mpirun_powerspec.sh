#!/bin/bash
# MOAB/Torque submission script for SciNet GPC 
#
#PBS -l nodes=2:ppn=8,walltime=1:00:00
#PBS -N powerspec
 
 
# DIRECTORY TO RUN - $PBS_O_WORKDIR is directory job was submitted from
cd $PBS_O_WORKDIR
 
# EXECUTION COMMAND; -np = nodes*ppn
mpirun -np 16 /home/r/rbond/gstein/src/powerspectra/powerspectra /home/r/rbond/gstein/src/powerspectra/data/6Gpc_n3072_nb16_nt16_13599_merge.pksc 6Gpc_2048 6000 2048