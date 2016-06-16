#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

rsync --progress -a -L ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_fastq_correct_merged/RIL*/*.gz ./

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

