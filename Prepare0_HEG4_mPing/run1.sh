#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
#PBS -V

cd $PBS_O_WORKDIR

bash /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4_FC70.1/run_these_jobs.sh

echo "Done"
