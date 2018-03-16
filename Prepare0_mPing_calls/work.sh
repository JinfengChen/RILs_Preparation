echo "RIL43 unique mPing"
grep "RIL43" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > RIL43.unique_mPing.gff
#22
bedtools window -w 10 -a RIL43.unique_mPing.gff -b RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.ALL.clean.gff | wc -l
#1 ping and 7 pong
less -S ~/BigData/00.RD/CNVRepeat_RIL/CNVRepeat_input_result_ping_random.summary1.txt
less -S ~/BigData/00.RD/CNVRepeat_RIL/CNVRepeat_input_result_pong_random.summary1.txt

echo "allele frquency for method" #20180315
grep "Parental" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency | awk '$3-$2 > 10'| cut -f7 | perl ~/BigData/software/bin/numberStat.pl

