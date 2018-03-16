ln -s ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.list
#10
grep "Parental" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency | awk '$7<0.3'| wc -l
#RIL list
grep -v "Sample" ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.txt | cut -f1 | sed 's/RIL//' > RILs.list
#shared_parental
grep "Parental" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency | awk '$7<0.3' | cut -f4 > Shared_Parental_0.3.list
#shared_ril
grep "RIL" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency | awk '$7>0.1' | cut -f4 > Shared_RIL_0.1.list

python heatmap_matrix.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.list
cat Shared_Parental_0.3.R| R --slave
cat Shared_RIL_0.1.R | R --slave

cat Shared_mPing.R| R --slave
