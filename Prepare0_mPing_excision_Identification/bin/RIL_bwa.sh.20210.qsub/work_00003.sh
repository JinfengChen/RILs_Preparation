perl /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Identification/bin/step1_Mapping.pl -ref /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Pseudogenome/MSU_r7.Pseudo_mPing_RILs.fa -1 /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Identification/input/RILs_ALL_unmapped_mping_fastq/RIL100/RIL100_1.fq -2 /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Identification/input/RILs_ALL_unmapped_mping_fastq/RIL100/RIL100_2.fq -min 0 --max 500 -cpu 12 --tool bwa --project /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Identification/bin/RIL100; echo This-Work-is-Completed!
