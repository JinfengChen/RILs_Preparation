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
#analysis mping copy number and read depth on ping vs. mping correlation
cut -f2 RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt | perl ~/BigData/software/bin/numberStat.pl
#mean=193, -10%=173 and +10%=212
python Sum_Ping_mPing_Narrow_range.py --code RIL275_RelocaTE.sofia.ping_code.table --output RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.narrow_range
awk '$2>=173 && $2<=212' RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.narrow_range.txt
python Sum_Ping_mPing_High_depth.py --code RIL275_RelocaTE.sofia.ping_code.table --output RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.high_depth
python Sum_Ping_mPing_High_depth_table.py --table RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt --output high_depth
python Sum_Ping_mPing_High_depth.py --code RIL275_RelocaTE.sofia.ping_code.table --output RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.high_narrow --narrow
python Sum_Ping_mPing_High_depth_table.py --table RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt --output high_narrow --narrow

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
ln -s ../Prepare0_mPing_excision_footprint/bin/Excision_newpipe_version1.footprint.list.noPing.txt ./
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency 
grep "Parental" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.parental.frequency 
grep "RIL" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.nonparental.frequency
echo "frequency figure"
cat mping.allele_frq.ALL_shared.R | R --slave
cat mping.allele_frq.parental.R | R --slave
cat mping.allele_frq.nonparental.R | R --slave
python mping.excision.avg_frequency.py
cat mping.excision.avg_frequency.R | R --slave


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
cd Prepare0_mPing_excision_Identification/bin
python runRIL_bwa.py --input ../input/RILs_ALL_unmapped_mping_fastq > log 2>&1 &
perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl -q js --maxjob 5 --lines 1 --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no RIL_bwa.sh > log1 2>&1 &
mkdir ../input/RILs_ALL_unmapped_mping_bam
mv *.bam ../input/RILs_ALL_unmapped_mping_bam/
mv *.bai ../input/RILs_ALL_unmapped_mping_bam/
mv *.dupli ../input/RILs_ALL_unmapped_mping_bam/
echo "3.1.4 summary mping boundary coverage"
cd Prepare0_mPing_excision_Identification/bin
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_275line_correct/MPR.geno.data ../input
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_275line_correct/MPR.geno.data ../input
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam_correct_merged --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam --gff_ref ../input/HEG4.ALL.mping.non-ref.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.gff > log 2>&1 &
#add ping code and excision/gt status to csv
python Ping_number_RILs.High_exicison.py --csv mPing_boundary_mPing --ping_code ../../Prepare0_mPing_calls/RIL275_RelocaTE.sofia.ping_code.table
#summary excision from csv
python Sum_excision_distance.py --dir mPing_boundary_mPing_GT_Ping_code --distance ../../Prepare0_mPing_distance/mPing_dist2.50Mb.list.sorted --blacklist Bam.Core.blacklist --project mPing_boundary.linked_50Mb_debug2 > log 2>&1 &
echo "3.2 Footprint"
mkdir Prepare0_mPing_excision_footprint
cd Prepare0_mPing_excision_footprint
mkdir bin
mkdir input
cd bin
python footprint_events.py --input ../../Prepare0_mPing_excision_Identification/bin/mPing_boundary.linked_50Mb_debug2.mping_excision.list --blacklist ../../Prepare0_mPing_excision_Identification/bin/Bam.Core.blacklist --output Excision_newpipe_version1 > log 2>&1 &
awk '$2>=5' Excision_newpipe_version1.footprint.list.txt > Excision_newpipe_version1.footprint.high.txt
cp Excision_newpipe_version1.footprint.list.txt Excision_newpipe_version1.footprint.list.noPing.txt

echo "3.3 Excision frequency plot"
cd Figure4_high_exicision_loci
ln -s ../Prepare0_mPing_excision_footprint/bin/Excision_newpipe_version1.footprint.list.txt ./
cut -f1-2 Excision_newpipe_version1.footprint.list.txt | grep "Chr" > Excision_newpipe_version1.footprint.list.draw.txt
cut -f2 Excision_newpipe_version1.footprint.list.draw.txt | perl ~/BigData/software/bin/numberStat.pl
cat mping.excision_events.binomial_test.R | R --slave
cat mping.excision_events.distr.R | R --slave

echo "3.4 mPing distance"
#need to finalize the distance, rewrite? How to get this "mPing_dist_RIL_AF0.1.50Mb.list.sorted" file? how many changes comapred to "mPing_dist2.50Mb.list.sorted"
#HEG4 distance
perl mPing_dist.pl --input HEG4.ALL.mping.non-ref.AF0.1.gff
sort -k3,3n -k1,1n -k2,2n mPing_dist2.50Mb.list | uniq > mPing_dist2.50Mb.list.sorted
#RIL distance
perl mPing_dist.pl --input RIL275_RelocaTEi.CombinedGFF.characterized.AF0.1.gff
sort -k3,3n -k1,1n -k2,2n mPing_dist_RIL_AF0.1.50Mb.list | uniq > mPing_dist_RIL_AF0.1.50Mb.list.sorted
#excision distance matrix
python Excision_Distance.py --excision1 Excision_newpipe_version1.footprint.list.txt --distance mPing_dist2.50Mb.list.sorted --gff HEG4.ALL.mping.non-ref.AF0.1.gff --output Excision_distance_HEG4.matrix_events
python Excision_Distance.py --excision1 Excision_newpipe_version1.footprint.list.txt --distance mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff HEG4.ALL.mping.non-ref.AF0.1.gff --output Excision_distance_RIL.matrix_events
paste Excision_distance.matrix_events.1.txt Excision_distance_HEG4.matrix_events.1.txt | awk -F"\t" '$3!=$6' | awk '$2>3'| less -S

echo "3.5 Excision and mPing distance"
cd Figure5_high_exicision_in_proximity
#excision vs. distance dotplot
ln -s ../Prepare0_mPing_distance/Excision_distance_HEG4.matrix_events.1.txt ./Excision_distance_HEG4.matrix_events.1.txt
ln -s ../Prepare0_mPing_distance/Excision_distance_RIL.matrix_events.1.txt ./Excision_distance_RIL.matrix_events.1.txt
cat Excision_distance.matrix_events.plot.R | R --slave
#high excision mPing vs. 413 parental mPing
python Excision_Distance_Boxplot.py --highexcision Excision_newpipe_version1.footprint.high.txt --distance mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff HEG4.ALL.mping.non-ref.AF0.1.gff --output Excision_Events_distance.boxplot
cat Excision_Events_distance.boxplot.R | R --slave

echo "3.6 Excision and mPing frequency"
mkdir Figure4_high_exicision_mping_frequency 
ln -s ../Prepare0_mPing_excision_footprint/bin/Excision_newpipe_version1.footprint.list.noPing.txt ./
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency ./
grep "Parental" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.parental.frequency
grep -v "Parental" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.nonparental.frequency

echo "3.7 Excision and Ping copy number"
mkdir Figure4_high_exicision_ping_cor 
python Excision_Number_Ping.py --excision Excision_newpipe_version1.footprint.list.noPing.rils.txt --ping RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.narrow_range.txt --output RIL272.RIL_mPing_Ping_Excision.narrow_range.table.txt
python Excision_Number_Ping.py --excision Excision_newpipe_version1.footprint.list.noPing.rils.txt --ping RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt --output RIL272.RIL_mPing_Ping_Excision.table.txt
cat ping_number_clean_excision_cor_plot.R | R --slave

