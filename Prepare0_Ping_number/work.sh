echo "sofia verison"
sed 's/RIL//' RIL275_RelocaTE.sofia.ping_code.table.sorted | sort -k3,3n| awk '{print "RIL"$3"\t"$2"\t"$1}' > RIL275_RelocaTE.sofia.ping_code.table.sorted.txt

echo "genotype"
python Ping_genotype.py --input Ping.list --snp_map RIL_275_correct/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_correct --output RILs_275_correct_Ping_genotype

echo "lulu"
#RIL_ping_copynumber_lulu_272

echo "RelocaTE2"
python Ping_RelocaTE2.py --input RILs_ALL_fastq_correct_merged_RelocaTEi_Ping > log 2>&1 &
sed 's/_RelocaTEi//' RILs_ALL_fastq_correct_merged_RelocaTEi_Ping.summary | sed 's/RIL//'| sort -k1,1n > RILs_ALL_fastq_correct_merged_RelocaTEi_Ping.summary.sorted.txt
sed 's/_RelocaTEi//' RILs_ALL_fastq_correct_merged_RelocaTEi_Ping.summary | sed 's/RIL//'| sort -k1,1n| grep "*"

echo "sofia vs lulu"
#65 diff
paste RIL_ping_copynumber_lulu_272.txt RIL275_RelocaTE.sofia.ping_code.table.sorted.txt | awk '$2 != $6' | wc -l
#12 might diff
paste RIL_ping_copynumber_lulu_272.txt RIL275_RelocaTE.sofia.ping_code.table.sorted.txt | awk '$2 == $6' | awk '$3!=$5' | wc -l

echo "genotype vs lulu"
#102 may diff
paste RIL_ping_copynumber_lulu_272.txt RILs_275_correct_Ping_genotype.table.ping_code.list | awk '$2!=$7'| wc -l
#2 may diff
paste RIL_ping_copynumber_lulu_272.txt RILs_275_correct_Ping_genotype.table.ping_code.list | awk '$2 == $7 && $3!=$8'| wc -l   
#check detail comfirmed 5 might be different

echo "new ping bam files"
python subbam_RILs.py --input New_Ping.list

