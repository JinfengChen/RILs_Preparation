#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=run_excision.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam_correct_merged --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing --gff_ref ../input/Parent.ALL.mPing.415.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing_415.gff --project mPing_boundary_415_nonref

echo "Done"
