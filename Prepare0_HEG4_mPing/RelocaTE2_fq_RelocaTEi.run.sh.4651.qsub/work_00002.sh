#!/bin/bash
python /rhome/cjinfeng/BigData/00.RD/RelocaTE2/scripts/relocaTE2.py --mate_1_id _p1 --mate_2_id _p2 --split --te_fasta /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa --genome_fasta /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU_r7.fa --fq_dir /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq/HEG4.2.2 --outdir /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi --reference_ins /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU_r7.fa.mPing.RepeatMasker.out --size 500 --step 1234567 --mismatch 2 --run --cpu 8 --aligner blat --verbose 3; echo This-Work-is-Completed!
