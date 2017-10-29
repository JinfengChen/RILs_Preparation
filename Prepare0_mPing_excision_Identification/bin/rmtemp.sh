#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=rmtemp.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

FILE=`ls *.sort.bam | head -n $N | tail -n 1`
sort_bam=$FILE
prefix=${FILE%.sort.bam}
bam=$prefix\.bam

echo "$bam"
echo "$sort_bam"

rmdup="/opt/linux/centos/7.x/x86_64/pkgs/picard/1.130/bin/picard MarkDuplicates"

rm $prefix\.1.bwa.log2 $prefix\.2.bwa.log2 $prefix\.1.sai $prefix\.2.sai $prefix\.convert.* $prefix\.raw.bam $prefix\.rmdup.* $prefix\.sam* $prefix\.sort.*
#$rmdup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT INPUT=$sort_bam OUTPUT=$bam METRICS_FILE=$prefix\.dupli

echo "Done"
