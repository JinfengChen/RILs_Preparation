echo "mPing calls from RelocaTE2"
mkdir Prepare0_mPing_calls
cd Prepare0_mPing_calls
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.* ./
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.summary ./
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.summary_clean.table ./
#generate clean gff of mPing call: remove TSD that not 3 bp and Ping/Pong calls.
python Clean_Calls.py --gff RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.gff
#generate unique gff of mPing call; generate unique mPing number according to different copy number of Ping
python Unique_mPing_clean.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.gff --reference HEG4.ALL.mping.non-ref.AF0.1.gff --code RIL275_RelocaTE.sofia.ping_code.table
#analyze shared mPing with HEG4 or in the RILs and unique mPing in each RILs
#generate summary table for shared and unique mPing in each RIL
python Sum_unique_clean.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.summary
#merge shared_unique_table with Ping code
python MergePingCode.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.txt
#generate *.class.summary which have summary of mPing in parental/shared/unique and hom/het/som
python Sum_class_clean.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff
#generate *.type.summary which have summary of hom/het/som for unique mPing
python Sum_type_clean.py --prefix RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean
#generate unique hom/het mPing number according to different copy number Ping
python Sum_Ping_mPing.py --code RIL275_RelocaTE.sofia.ping_code.table --output RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean
echo "Fig1a"
cd Figure1_mPing_insertions
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.class.summary ./
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.type.summary ./
cat Fig1a.R | R --slave
echo "Fig1b"
cd Figure2_mPing_Ping_correlation
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.ping_number.summary ./
cat Fig1b.R | R --slave

echo "mPing frequency"
cd Figure3_mPing_frequency
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency 
grep "Parental" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency 
grep "RIL" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.nonparental.frequency
echo "frequency figure"
cat mping.allele_frq.ALL_shared.R | R --slave
cat mping.allele_frq.parental.R | R --slave
cat mping.allele_frq.nonparental.R | R --slave

echo "3 mPing excision"
echo "3.1 Identification"
mkdir Prepare0_mPing_excision_Identification
cd Prepare0_mPing_excision_Identification
mkdir bin
mkdir input
cd input
ln -s ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_correct_merged/ ./
#Parent.ALL.mPing.gff includes 51 reference mPing, 7 HEG4 Ping, 1 NB Ping and 545 HEG4 mPing.
cp ~/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/Parent.ALL.mPing.gff ./ 
echo "3.1.1 Prepare unmapped and mPing region reads"
cd bin
python PrepareFastq.py --bam ../input/RILs_ALL_bam_correct_merged > log 2>&1 &
echo "3.1.2 Prepare pseduogenome with mPing/ping inserted"
mkdir Prepare0_mPing_excision_Pseudogenome
cd Prepare0_mPing_excision_Pseudogenome 
python Pseudo_TEinsertion_Genome.py --repeat mPing_Ping_Pong.fa --gff HEG4.ALL.mping.non-ref.gff --genome MSU_r7.fa
python Pseudo_TEinsertion_Genome.py --repeat mPing_Ping_Pong.fa --gff RIL275_RelocaTEi.CombinedGFF.characterized.AF0.1.gff --genome MSU_r7.fa --project MSU_r7.Pseudo_mPing_RILs
/opt/linux/centos/7.x/x86_64/pkgs/bwa/0.7.12/bin/bwa index MSU_r7.Pseudo_mPing.fa
/opt/linux/centos/7.x/x86_64/pkgs/bwa/0.7.12/bin/bwa index MSU_r7.Pseudo_mPing_RILs.fa
echo "3.1.3 Mapping reads to pseudogenome"
python runRIL_bwa.py --input ../input/RILs_ALL_unmapped_mping_fastq > log 2>&1 &
perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl -q js --maxjob 5 --lines 1 --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no RIL_bwa.sh > log1 2>&1 &
mkdir ../input/RILs_ALL_unmapped_mping_bam
mv *.bam ../input/RILs_ALL_unmapped_mping_bam/
mv *.bai ../input/RILs_ALL_unmapped_mping_bam/
mv *.dupli ../input/RILs_ALL_unmapped_mping_bam/

