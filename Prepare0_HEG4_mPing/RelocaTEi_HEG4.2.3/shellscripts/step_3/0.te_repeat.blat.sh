/usr/local/bin/blat -minScore=10 -tileSize=7 /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping.fa /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA_1.fa /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/blat_output/HEG4.2.3.MSU7_BWA_1.te_repeat.blatout 1>>/rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/blat_output/blat.out 2>>/rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/blat_output/blat.out
python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE_trim.py /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/blat_output/HEG4.2.3.MSU7_BWA_1.te_repeat.blatout /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA_1.fq 10 1 > /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/flanking_seq/HEG4.2.3.MSU7_BWA_1.te_repeat.flankingReads.fq
