cat /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/results/*.all_nonref_insert.gff > /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/results/ALL.all_nonref_insert.gff
cat /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/results/*.all_nonref_insert.txt | grep "^TE" -v > /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/results/ALL.all_nonref_insert.txt
perl /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/characterizer.pl -s /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/results/ALL.all_nonref_insert.txt -b /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/bam/HEG4.2.3.MSU7_BWA.bam -g /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa --samtools /opt/samtools-0.1.16/samtools
