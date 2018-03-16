/opt/samtools-0.1.16/samtools view -h /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/bam/HEG4.2.3.MSU7_BWA.bam | awk '$5<60' | samtools view -Shb - | samtools sort -m 500000000 -n - /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA.subset 2> /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/run.std
/opt/bedtools/2.17.0-25-g7b42b3b/bin//bedtools bamtofastq -i /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA.subset.bam -fq /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA_1.fq -fq2 /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA_2.fq 2> /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/run.std
/rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE_fq2fa.pl /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA_1.fq /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA_1.fa
/rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE_fq2fa.pl /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA_2.fq /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3/repeat/fastq/HEG4.2.3.MSU7_BWA_2.fa
