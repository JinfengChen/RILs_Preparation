#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=Target_Run.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

module load ncbi-blast/2.2.26
module load mafft
module load treebest

export PATHONPATH=$PATHONPATH:/rhome/cjinfeng/BigData/software/pythonlib/lib/python2.7/site-packages:/rhome/cjinfeng/BigData/software/Target/2.0/
export PATH=$PATH:/rhome/cjinfeng/BigData/software/fasttree/
CPU1=$SLURM_NTASKS
CPU2=1
mode=mi

/rhome/cjinfeng/BigData/software/Target/2.0/target.py -q /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_Ping_variants/query/ping.fa -t nucl -o /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_Ping_variants/Target -i $mode -P $CPU1 -C $CPU2 -b_a 12000 -b_d 10 -E -W 5 -f 1.2 -a flanks -p_M 0.40 -p_n 12000 -p_d 6000 -p_f 300 -S 'MSA' /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_Ping_variants/HEG4_ALLPATHLG_v1.chr.fasta Target_Run_HEG4_ping


echo "Done"
