#!/bin/bash

# MOAB/Torque submission script for SciNet GPC   
#     
#PBS -l nodes=5:ppn=8,walltime=4:00:00     
#PBS -N powerspec             

# DIRECTORY TO RUN - $PBS_O_WORKDIR is directory job was submitted from       
cd $PBS_O_WORKDIR

nxi=150
folder=/scratch2/r/rbond/gstein/peak-patch-runs/current/PAPER_RUNS/FINAL_RUNS/
Lbox=1750
nrespk=2000 
fmt=0 # 0 for peaks, 1 for field
nproc=40

#OUTDATA
dirnamenp=numpy_data
if [ ! -d $dirnamenp ]; then mkdir $dirnamenp ; fi

#FIELD
dirnameF=ncutF_$Lbox_$nrespk
if [ ! -d $dirnameF ]; then mkdir $dirnameF; fi

python pk_comparison.py 1832.46 $nrespk 1    $nxi 0 0 0 $folder $dirnameF $nproc 
python data_append_pkxi.py $dirnameF 0

#HALO CATALOGUES
for ncut in `echo 10000 30000 100000 300000`
do
#    dirname=ncut_$Lbox_$nrespk_$ncut
    dirnameL=ncutL_$Lbox_$nrespk_$ncut

#    if [ ! -d $dirname  ]; then mkdir $dirname ; fi
    if [ ! -d $dirnameL ]; then mkdir $dirnameL; fi

#    python pk_comparison.py $Lbox $nrespk $fmt $nxi $ncut 0 0 $folder $dirname $nproc 
    python pk_comparison.py $Lbox $nrespk $fmt $nxi $ncut 1 0 $folder $dirnameL $nproc 

#    python data_append_pkxi.py $dirname $ncut
    python data_append_pkxi.py $dirnameL $ncut
done

