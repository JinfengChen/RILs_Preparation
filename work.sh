echo "mPing calls from RelocaTE2"
mkdir Prepare0_mPing_calls
cd Prepare0_mPing_calls
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.* ./
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.summary ./
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.summary_clean.table ./
#generate shared mPing wiht NB and both NB and HEG4 for each RIL
cp ~/BigData/00.RD/RILs/Transpostion/bin/Compare_NB_mPing/HEG4.mping.ref_only.gff ./
cp ~/BigData/00.RD/RILs/Transpostion/bin/Compare_NB_mPing/HEG4.mping.shared.gff ./
cp /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_Parental_Ping/HEG4.ALL_Filter.ping.gff ./
python Sum_Shared_mPing_Ref.py --gff RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.gff

#generate clean gff of mPing call: remove TSD that not 3 bp and Ping/Pong calls.
python Clean_Calls.py --gff RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.gff
#generate unique gff of mPing call; generate unique mPing number according to different copy number of Ping
#cp ../Prepare0_HEG4_mPing/HEG4.P.RelocaTE2.mping.non-ref.AF0.1.gff ./
#python Unique_mPing_clean.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.gff --reference HEG4.ALL.mping.non-ref.AF0.1.gff --code RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt
#python Unique_mPing_clean.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.gff --reference HEG4.P.RelocaTE2.mping.non-ref.AF0.1.gff --code RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt
cp ../Prepare0_HEG4_mPing/HEG4.2.3.RelocaTE2.mping.non-ref.AF0.1.gff
python Unique_mPing_clean.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.gff --reference HEG4.2.3.RelocaTE2.mping.non-ref.AF0.1.gff --code RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt
#use HEG4 and NB mPing as reference mPing: HEG4_NB_ALL.mPing.gff
#python Unique_mPing_clean.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.ALL_and_Shared.gff --reference HEG4_NB_ALL.mPing.gff --code RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt
#NB, cat HEG4.mping.ref_only.gff HEG4.mping.shared.gff > HEG4.mping.466.gff
python Unique_mPing_clean_NB.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.clean.gff --reference HEG4.mping.466.gff
ln -s RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.clean.overlap_ref_NB RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.overlap_ref_NB
ln -s RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.clean.overlap_ref_ril RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.overlap_ref_ril

#numbers
#non-reference mPing: 80109 mpings or 16865 mping loci 
wc -l RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.ALL.clean.gff
awk '{print $1":"$4"-"$5}' RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.ALL.clean.gff | sort | uniq | wc -l
#reference or shared mPing: 7345 mpings or 51 mping loci
grep -v "23521641" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.gff > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.clean.gff
wc -l RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.clean.gff
awk '{print $1":"$4"-"$5}' RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.clean.gff | sort | uniq | wc -l
#total: 87454 mping or 16916 mping loci
cat RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.ALL.clean.gff RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.clean.gff > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.ALL_and_Shared.clean.gff
awk '{print $1"-"$4"-"$5}' RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.ALL_and_Shared.clean.gff| sort | uniq | wc -l
ln -s RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.ALL_and_Shared.clean.gff RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.ALL_and_Shared.gff
#HEG4_NB all mPing: 466=415+51
cat HEG4.2.3.RelocaTE2.mping.non-ref.AF0.1.gff HEG4.mping.ref_only.gff HEG4.mping.shared.gff > HEG4_NB_ALL.mPing.gff

#analyze shared mPing with HEG4 or in the RILs and unique mPing in each RILs
#generate summary table for shared and unique mPing in each RIL
#generate mPing frequency table
python Sum_unique_clean.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.summary
#merge shared_unique_table with Ping code
python MergePingCode.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.txt
#generate *.class.summary which have summary of mPing in parental/shared/unique and hom/het/som
ln -s RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.ALL_and_Shared.clean.gff RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.ALL_and_Shared.gff
python Sum_class_clean.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff

#12675 mping or 1915 mping loci
awk '{print $1":"$4"-"$5}' RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_nonparental.gff | sort | uniq | wc -l
#60244 mping or 466 mping loci
awk '{print $1":"$4"-"$5}' RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_parental.gff | sort | uniq | wc -l
#14535 uniq mping
wc -l RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff

#generate *.type.summary which have summary of hom/het/som for unique mPing
python Sum_type_clean.py --prefix RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean
#generate unique hom/het mPing number according to different copy number Ping
#python Sum_Ping_mPing.py --code RIL275_RelocaTE.sofia.ping_code.table --output RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean
python Sum_Ping_mPing.py --code RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt --output RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean
#analysis mping copy number and read depth on ping vs. mping correlation
cut -f2 RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt | perl ~/BigData/software/bin/numberStat.pl
#mean=193, -10%=173 and +10%=212
python Sum_Ping_mPing_Narrow_range.py --code RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt --output RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.narrow_range
awk '$2>=173 && $2<=212' RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.narrow_range.txt
python Sum_Ping_mPing_High_depth.py --code RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt --output RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.high_depth
python Sum_Ping_mPing_High_depth_table.py --table RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt --output high_depth
python Sum_Ping_mPing_High_depth.py --code RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt --output RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.high_narrow --narrow
python Sum_Ping_mPing_High_depth_table.py --table RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt --output high_narrow --narrow
#add ref_only and shared mPing
paste RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.Shared.ref_shared.table.txt | cut -f1-10,12- > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.txt
grep -v "Sample" RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.txt | awk '{print $2+$11}' | perl ~/BigData/software/bin/numberStat.pl
#mean=221, -10%=199 and +10%=243 
awk '$2>=199 && $2<=243' RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.txt > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.narrow_range.txt
python Sum_Ping_mPing_High_depth_table.py --table RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.txt --output high_depth
python Sum_Ping_mPing_High_depth_table.py --table RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.txt --output high_narrow --narrow


echo "Fig1a"
cd Figure1_mPing_insertions
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.class.summary ./
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.type.summary ./
cat Fig1a.R | R --slave
echo "Fig1b"
cd Figure2_mPing_Ping_correlation
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.ping_number.summary ./
cat Fig1b.R | R --slave
#single ping analysis
paste RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt ../Prepare0_Sequence_Depth/RILs_ALL_bam_correct_merged.sorted.summary > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.depth.txt
head -n 1 RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.depth.txt > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.depth.single_ping.txt
awk '$9==1' RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.depth.txt >> RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.depth.single_ping.txt

echo "Fig2"
#Fig2c
cd Figure0_QTL
mkdir QTL_unique_mPing
cd QTL_unique_mPing
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_272line/MPR.cross.uniq.cro ./
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_272line/MPR.cross.uniq.map ./
cat MPR.cross.uniq_mPing_plot.R| R --slave
#Plot Bin and Tree
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_272line/MPR.geno.bin.uniq ./
cp ~/BigData/00.RD/RILs/Figures/QTL_BinPlot/bin/MapDist.py ./
cp ~/BigData/00.RD/RILs/Figures/QTL_BinPlot/bin/MapGenome_Chr.py
cp ~/BigData/00.RD/RILs/Figures/QTL_BinPlot/bin/QTL_BinMap_Chr.py ./
cp ~/BigData/00.RD/RILs/Figures/QTL_BinPlot/bin/TraitDist4Tree.py ./
cp ~/BigData/00.RD/RILs/QTL_pipe/input/trait/May28_2013.RIL.trait.table.QTL.trait.txt.first ./
python MapDist.py --input MPR.geno.bin.uniq > MPR.geno.bin.uniq.dist
python MapGenome_Chr.py --input MPR.geno.bin.uniq > MSU7.Chr.midpoint
python QTL_BinMap_Chr.py --input MPR.geno.bin.uniq.dist --bin MPR.geno.bin.uniq.new
python TraitDist4Tree.py --input May28_2013.RIL.trait.table.QTL.trait.txt.first > May28_2013.RIL.trait.table.QTL.trait.txt.first.dist
python Trait_Tree_Label.py --input May28_2013.RIL.trait.table.QTL.trait.txt.first.dist
mv Trait_Tree_Label.pdf Trait_Tree_Label.All_Traits.pdf
python Trait_Tree_Label.py --input MPR.geno.bin.uniq.dist
mv Trait_Tree_Label.pdf Trait_Tree_Label.Bin_Map.pdf
#summary map
mkdir MapSummary
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_272line/MPR.geno.bin.uniq ./
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_272line/MPR.cross.uniq.cro ./
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_272line/MPR.cross.uniq.map ./
python MapSummary.py --input MPR.cross.uniq --bin MPR.geno.bin.uniq > RIL272.log

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
#20180302, Pseudogenome for 51 Reference and Shared mPing in HEG4 and NB
python Pseudo_TEinsertion_Genome_Ref.py --gff Parent.ALL.mPing.Ref_Shared.gff --genome MSU_r7.fa --project Parent.Pseudo_mPing.Ref_Shared

echo "3.1.3 Mapping reads to pseudogenome"
cd Prepare0_mPing_excision_Identification/bin
python runRIL_bwa.py --input ../input/RILs_ALL_unmapped_mping_fastq > log 2>&1 &
perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl -q js --maxjob 5 --lines 1 --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no RIL_bwa.sh > log1 2>&1 &
mkdir ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing
mkdir ../input/RILs_ALL_unmapped_mping_bam_RILs_mPing
#run and mv using HEG4 mPing
mv *.bam ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing/
mv *.bai ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing/
mv *.dupli ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing/
#run and mv using RILs mPing > AF0.1
mv *.bam ../input/RILs_ALL_unmapped_mping_bam_RILs_mPing/
mv *.bai ../input/RILs_ALL_unmapped_mping_bam_RILs_mPing/
mv *.dupli ../input/RILs_ALL_unmapped_mping_bam_RILs_mPing/
#20180303, mapped to Pseudogenome
python runRIL_bwa.py --input ../input/RILs_ALL_unmapped_mping_fastq > log 2>&1 &
#perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 30 --lines 1 --interval 120 --task 12 --mem 20G --time 100:00:00 --convert no RIL_bwa.sh > log1 2>&1 &
perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 10 --queue intel --lines 1 --interval 120 --task 12 --mem 20G --time 100:00:00 --convert no RIL_bwa_1_100.sh > log1 2>&1 &
perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 10 --queue batch  --lines 1 --interval 120 --task 12 --mem 20G --time 100:00:00 --convert no RIL_bwa_101_199.sh > log2 2>&1 &
perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 10 --queue stajichlab --lines 1 --interval 120 --task 12 --mem 20G --time 100:00:00 --convert no RIL_bwa_200_272.sh > log3 2>&1 &

mkdir ../input/RILs_ALL_unmapped_mping_bam_Ref_mPing
mv *.bam ../input/RILs_ALL_unmapped_mping_bam_Ref_mPing/
mv *.bai ../input/RILs_ALL_unmapped_mping_bam_Ref_mPing/
mv *.dupli ../input/RILs_ALL_unmapped_mping_bam_Ref_mPing/

echo "3.1.4 summary mping boundary coverage"
cd Prepare0_mPing_excision_Identification/bin
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_275line_correct/MPR.geno.data ../input
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_275line_correct/MPR.geno.data ../input
#RILs mPing, used to design PCR screen
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam_correct_merged --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam_RILs_mPing --gff_ref ../input/RIL275_RelocaTEi.CombinedGFF.characterized.AF0.1.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing_RILs.gff > log 2>&1 &
python Ping_number_RILs.High_exicison.py --csv mPing_boundary_mPing_RILs --ping_code ../../Prepare0_mPing_calls/RIL275_RelocaTE.sofia.ping_code.table

#HEG4 mPing, used to call excision
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam_correct_merged --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing --gff_ref ../input/HEG4.ALL.mping.non-ref.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.gff > log 2>&1 &
#add ping code and excision/gt status to csv
python Ping_number_RILs.High_exicison.py --csv mPing_boundary_mPing --ping_code ../../Prepare0_mPing_calls/RIL275_RelocaTE.sofia.ping_code.table
#summary excision from csv
python Sum_excision_distance.py --dir mPing_boundary_mPing_GT_Ping_code --distance ../../Prepare0_mPing_distance/mPing_dist2.50Mb.list.sorted --blacklist Bam.Core.blacklist --project mPing_boundary.linked_50Mb_debug2 > log 2>&1 &

#51 Ref and Shared mPing, 20180316
python mPing_Boundary_Coverage_Ref.py --bam_ref ../input/RILs_ALL_unmapped_mping_bam_Ref_mPing/ --bam_pseudo ../input/RILs_ALL_bam_correct_merged/ --gff_ref ../input/Parent.Pseudo_mPing.Ref_Shared.gff --gff_pseudo ../input/Parent.ALL.mPing.Ref_Shared.gff > log 2>&1 &
python Ping_number_RILs.High_exicison_Ref.py --csv mPing_boundary_mPing --ping_code ../../Prepare0_mPing_calls/RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt

#merge 415 non_ref and 51 ref_shared
cp mPing_boundary_51_ref_mPing_GT_Ping_code mPing_boundary_mPing_GT_Ping_code
cp -R mPing_boundary_415_nonref_mPing_GT_Ping_code mPing_boundary_mPing_GT_Ping_code
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
echo "summary excision number"
python sum_excision.py --input Excision_newpipe_version1.footprint.list.noPing.txt > Excision_newpipe_version1.footprint.list.noPing.summary
python sum_excision_Table_S.py --input Excision_newpipe_version1.footprint.list.noPing.txt

echo "3.3 Excision frequency plot"
cd Figure4_high_exicision_loci
ln -s ../Prepare0_mPing_excision_footprint/bin/Excision_newpipe_version1.footprint.list.noPing.txt ./
cut -f1-2 Excision_newpipe_version1.footprint.list.noPing.txt | grep "Chr" > Excision_newpipe_version1.footprint.list.draw.txt
cut -f2 Excision_newpipe_version1.footprint.list.draw.txt | perl ~/BigData/software/bin/numberStat.pl
cat mping.excision_events.binomial_test.R | R --slave
cat mping.excision_events.distr.R | R --slave
cd Prepare0_Parental_mPing_loci
python Mendel_Deviation.py --input RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency | sort -k9nr | awk '$9<0.0001' | sort -k7n > RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.deviation_from_mendel.txt

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
#high excision distance
python High_excision_distance.py --input 1 > Excision_newpipe_version1.footprint.high.distance.txt

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
#add Ref_only and shared mPing
python Excision_Number_Ping_AddRefmPing.py --excision Excision_newpipe_version1.footprint.list.noPing.rils.txt --ping RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.narrow_range.txt --output RIL272.RIL_mPing_Ping_Excision.Ref_mPing.narrow_range.table.txt
python Excision_Number_Ping_AddRefmPing.py --excision Excision_newpipe_version1.footprint.list.noPing.rils.txt --ping RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.Ref_mPing.txt --output RIL272.RIL_mPing_Ping_Excision.Ref_mPing.table.txt
cat ping_number_clean_excision_cor_plot_AddRefmPing.R | R --slave

echo "3.8" High excision associated genes
mkdir Prepare0_mPing_excision_genes
cd Prepare0_mPing_excision_genes 
python High_excision_asso_genes.py --input 1 > Excision_newpipe_version1.footprint.high.associated_gene.txt

echo "3.9" Excision Validation
echo "3.9.1"
mkdir Prepare0_mPing_SV_screen/bin
python mPing_unknown_gt_pairs_RILs.py --input High_excision_mPing.pairs --matrix ../input/mPing_boundary_mPing_RILs_GT_Ping_code
python RILs_need.py > High_excision_mPing.needDNA

echo "3.9.2" PCR screen for SV
#update blackout RILs that need for PCR.
mkdir Prepare0_mPing_SV_screen_update/bin
cd Prepare0_mPing_SV_screen_update/bin
cp ../../Prepare0_mPing_SV_screen/bin/*.check.list ./
cp ../../Prepare0_mPing_SV_screen/bin/*.rils ./
python RILs_need_from_blackout.py > ALL.blackout.update

#after pcr screen, summary the results
mkdir Prepare0_mPing_SV_screen_PCR/bin
cd Prepare0_mPing_SV_screen_PCR/bin 
python SV_screen_PCR_results.py
python SV_screen_PCR_not_finished.py > RILs.not_finished
grep -v "high" RILs.not_finished | grep -v "HEMP" | cut -f1 | sort | uniq | sort -n > RILs.not_finished.DNA

echo "4. mPing distribution"
mkdir Figure7_Circos 
ln -s ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff ./
grep "hom" ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > mPing_gff/RIL.gff
grep -e "het" -e "som" ../Prepare0_mPing_calls/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > mPing_gff/Somatic.gff
cp ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Figures/Fig2_Distribution/mPing_gff/Strains.gff mPing_gff/
cp ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Figures/Fig2_Distribution/mPing_gff/Simulate0001.gff mPing_gff/
python /rhome/cjinfeng/BigData/00.RD/Circos/bin/scripts/distri_data_pre_gff.py --head MSU7.circos.head --input mPing.gff.RelocaTE2.list
cp ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Figures/Fig2_Distribution/mPing_Distribution_Circos/OpenChromatin.gff.list ./
python /rhome/cjinfeng/BigData/00.RD/Circos/bin/scripts/distri_data_pre_gff.py --head MSU7.circos.head --input OpenChromatin.gff.list

#preparing mping distribution around gene
mkdir Figure7_mPing_distribution_pre Figure7_mPing_distribution_pre/bin Figure7_mPing_distribution_pre/input
ln -s ../../Figure7_Circos/mPing_gff/ ./
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Figures/Fig2_Distribution/mPing_gff/MSU_r7.all.final.* mPing_gff/ 
#Somatic
bedtools intersect -a ./mPing_gff/Somatic.gff -b ./mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > Somatic.intersect
bedtools closest -a ./mPing_gff/Somatic.gff -b ./mPing_gff/MSU_r7.all.final.mRNA.gff -d > Somatic.mRNA.intersect
python mPing_position.py --input Somatic.intersect --mrna Somatic.mRNA.intersect
python mPing_intron.py --input Somatic.intersect
python mPing_intergenic.py --input Somatic.mRNA.intersect
python mPing_position_2kb.py --input Somatic.intersect --mrna Somatic.mRNA.intersect
#RIL
bedtools intersect -a ./mPing_gff/RIL.gff -b ./mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > RIL.intersect
bedtools closest -a ./mPing_gff/RIL.gff -b ./mPing_gff/MSU_r7.all.final.mRNA.gff -d > RIL.mRNA.intersect
python mPing_position.py --input RIL.intersect --mrna RIL.mRNA.intersect
python mPing_position_2kb.py --input RIL.intersect --mrna RIL.mRNA.intersect
python mPing_intron.py --input RIL.intersect
python mPing_intergenic.py --input RIL.mRNA.intersect
#Strain
bedtools intersect -a ./mPing_gff/Strains.gff -b ./mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > Strains.intersect
bedtools closest -a ./mPing_gff/Strains.gff -b ./mPing_gff/MSU_r7.all.final.mRNA.gff -d > Strains.mRNA.intersect
python mPing_position.py --input Strains.intersect --mrna Strains.mRNA.intersect
python mPing_position_2kb.py --input Strains.intersect --mrna Strains.mRNA.intersect
python mPing_intron.py --input Strains.intersect 
python mPing_intergenic.py --input Strains.mRNA.intersect 
#Simulation
bedtools intersect -a ./mPing_gff/Simulate0001.gff -b ./mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > Simulation.intersect
bedtools closest -a ./mPing_gff/Simulate0001.gff -b ./mPing_gff/MSU_r7.all.final.mRNA.gff -d > Simulation.mRNA.intersect
python mPing_position.py --input Simulation.intersect --mrna Simulation.mRNA.intersect
python mPing_position_2kb.py --input Simulation.intersect --mrna Simulation.mRNA.intersect
python mPing_intron.py --input Simulation.intersect 
python mPing_intergenic.py --input Simulation.mRNA.intersect


#chromome4 distribution
mkdir Figure7_mPing_distribution Figure7_mPing_distribution/bin Figure7_mPing_distribution/input Figure7_mPing_distribution/input/Chromosome4 
grep "Chr4" ../../Figure7_Circos/MSU7.exon.histogram.txt > ../input/Chromosome4/Chr4.exon.histogram.txt
grep "Chr4" ../../Figure7_Circos/GFF.Simulation.histogram.txt > ../input/Chromosome4/Chr4.simulation.histogram.txt
grep "Chr4" ../../Figure7_Circos/GFF.RIL.histogram.txt > ../input/Chromosome4/Chr4.RIL.histogram.txt
grep "Chr4" ../../Figure7_Circos/GFF.DHS.histogram.txt > ../input/Chromosome4/Chr4.DHS.histogram.txt
cat Plot_Density_Chr4.R | R --slave
#mping distance to DHS plot
python Distance2DHS.py > log 2>&1 &
#DHS profile around mping TSD
python TSDprofile_DHS_R.py --bam ../input/DHS.unique.bam --gff ../input/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > log 2>&1 &
#mPing distribution around gene
cat mping_intergenic_3distance_withsim_chromatin.R | R --slave
cat mping_intergenic_5distance_withsim_chromatin.R | R --slave
cat mping_intergenic_5distance_withsim_chromatin.R | R --slave
cat mping_position_breakY_withsim_chromatin.R | R --slave
cat mping_position_breakY_withsim_chromatin_2kb.R | R --slave
#DHS and Nucleosome distribution at TSS and TTS
qsub run_DHS.sh
qsub run_Nucleosome.sh
cat TSSprofile_Nuc_DHS.R | R --slave
cat TTSprofile_Nuc_DHS.R | R --slave

