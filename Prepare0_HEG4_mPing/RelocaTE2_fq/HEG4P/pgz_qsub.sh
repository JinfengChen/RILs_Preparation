#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=20G
#SBATCH --time=10:00:00
#SBATCH --output=pgz.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

FILE=`ls *.fq | head -n $N | tail -n 1`
if [ ! -e $FILE\.gz ]; then
   pigz $FILE -p 16
   chmod 664 $FILE\.gz
fi

echo "Done"
