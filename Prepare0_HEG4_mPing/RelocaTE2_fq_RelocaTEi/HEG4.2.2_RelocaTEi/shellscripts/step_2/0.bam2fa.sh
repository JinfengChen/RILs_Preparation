/rhome/cjinfeng/BigData/00.RD/RelocaTE2/bin/samtools sort -m 1000000000 -n /rhome/cjinfeng/BigData/01.Rice_genomes/HEG4/00.Bam/HEG4_MSU7_BWA/HEG4_2.2.MSU7_BWA.bam /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/repeat/fastq/HEG4_2.2.MSU7_BWA.sortbyname 2> /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/run.std
/rhome/cjinfeng/BigData/00.RD/RelocaTE2/bin/bedtools bamtofastq -i /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/repeat/fastq/HEG4_2.2.MSU7_BWA.sortbyname.bam -fq /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/repeat/fastq/HEG4_2.2.MSU7_BWA_1.fq -fq2 /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/repeat/fastq/HEG4_2.2.MSU7_BWA_2.fq 2> /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/run.std
/rhome/cjinfeng/BigData/00.RD/RelocaTE2/bin/seqtk seq -A /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/repeat/fastq/HEG4_2.2.MSU7_BWA_1.fq > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/repeat/fastq/HEG4_2.2.MSU7_BWA_1.fa
/rhome/cjinfeng/BigData/00.RD/RelocaTE2/bin/seqtk seq -A /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/repeat/fastq/HEG4_2.2.MSU7_BWA_2.fq > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq_RelocaTEi/HEG4.2.2_RelocaTEi/repeat/fastq/HEG4_2.2.MSU7_BWA_2.fa
