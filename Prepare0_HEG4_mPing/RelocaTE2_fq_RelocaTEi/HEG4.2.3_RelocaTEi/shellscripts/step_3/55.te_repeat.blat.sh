/rhome/cjinfeng/BigData/00.RD/RelocaTE2/bin/blat -minScore=10 -tileSize=7 /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.3_RelocaTEi/repeat/fastq_split/p115.SRR1619147_p2.fa /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.3_RelocaTEi/repeat/blat_output/p115.SRR1619147_p2.te_repeat.blatout 1>>/bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.3_RelocaTEi/repeat/blat_output/blat.out 2>>/bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.3_RelocaTEi/repeat/blat_output/blat.out
python /rhome/cjinfeng/BigData/00.RD/RelocaTE2/scripts/relocaTE_trim.py /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.3_RelocaTEi/repeat/blat_output/p115.SRR1619147_p2.te_repeat.blatout /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.3_RelocaTEi/repeat/fastq_split/p115.SRR1619147_p2.fq 10 10 2 > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.3_RelocaTEi/repeat/flanking_seq/p115.SRR1619147_p2.te_repeat.flankingReads.fq
