#!/bin/bash

# MOAB/Torque submission script for SciNet GPC   
#     
#PBS -l nodes=1:ppn=8,walltime=4:00:00     
#PBS -N powerspec             

# DIRECTORY TO RUN - $PBS_O_WORKDIR is directory job was submitted from       
cd $PBS_O_WORKDIR

nxi=95
folder=/scratch2/r/rbond/gstein/peak-patch-runs/current/PAPER_RUNS/FINAL_RUNS/output/
seed=13579

Lbox=1750
nrespk=512 
fmt=0 # 0 for peaks, 1 for field
nproc=4

dirnamenp=numpy_data
if [ ! -d $dirnamenp]; then mkdir $dirnamenp ; fi

for ncut in `echo 10000 30000 100000 300000`
do

    dirnameL=ncutL_$Lbox_$nrespk_$ncut
    dirname=ncut_$Lbox_$nrespk_$ncut

    if [ ! -d $dirnameL ]; then mkdir $dirnameL; fi
    if [ ! -d $dirname  ]; then mkdir $dirname ; fi


 
    cd $dirname
    python ../pk_comparison.py $Lbox $nrespk $fmt $nxi $ncut 0 0 $folder $nproc &
    cd ../

    cd $dirnameL
    python ../pk_comparison.py $Lbox $nrespk $fmt $nxi $ncut 1 0 $folder $nproc 
    cd ../

    python data_append_pkxi.py $dirname $ncut
    python data_append_pkxi.py $dirnameL $ncut
done



