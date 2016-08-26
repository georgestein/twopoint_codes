#!/bin/bash

nxi=10
folder=/scratch2/r/rbond/gstein/peak-patch-runs/current/PAPER_RUNS/FINAL_RUNS/output/
seed=13579

Lbox=1750
nrespk=512 
fmt=0 # 0 for peaks, 1 for field
nproc=4
for ncut in `echo 10000`
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

done



