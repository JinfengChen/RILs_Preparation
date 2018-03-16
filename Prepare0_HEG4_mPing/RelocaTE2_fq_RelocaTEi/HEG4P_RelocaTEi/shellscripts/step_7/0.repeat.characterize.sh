cat /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/*.all_nonref_insert.gff > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/ALL.all_nonref_insert.gff
cat /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/*.all_nonref_insert.txt | grep "^TE" -v > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/ALL.all_nonref_insert.txt
python /rhome/cjinfeng/BigData/00.RD/RelocaTE2/scripts/clean_false_positive.py --input /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/ALL.all_nonref_insert.gff --refte /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/existingTE.bed --bedtools /rhome/cjinfeng/BigData/00.RD/RelocaTE2/bin/bedtools
cat /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/*.all_ref_insert.txt > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/ALL.all_ref_insert.txt
cat /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/*.all_ref_insert.gff > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/ALL.all_ref_insert.gff
perl /rhome/cjinfeng/BigData/00.RD/RelocaTE2/scripts/characterizer.pl -s /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/ALL.all_nonref_insert.txt -b /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/HEG4_P.bam -g /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU_r7.fa --samtools /rhome/cjinfeng/BigData/00.RD/RelocaTE2/bin/samtools
