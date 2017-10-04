#python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE.py --te_fasta /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping.fa --genome_fasta /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa --bam /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/bam/HEG4.2.3.MSU7_BWA.bam --outdir /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.3 --reference_ins /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa.RepeatMasker.out --size 180
#python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE.py --te_fasta /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping.fa --genome_fasta /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa --bam /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/bam/FC70.1.MSU7_BWA.bam --outdir /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4_FC70.1 --reference_ins /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa.RepeatMasker.out --size 400
#python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE.py --te_fasta /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping.fa --genome_fasta /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa --bam /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/bam/HEG4_P.bam --outdir /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4_P --reference_ins /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa.RepeatMasker.out --size 180


#bash /rhome/cjinfeng/software/tools/TEMP/scripts/TEMP_Insertion_qsub.sh -i /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/bam/HEG4.2.3.MSU7_BWA.bam -s /rhome/cjinfeng/software/tools/TEMP/scripts -r /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping_pogo.fa -t /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa.RepeatMasker.out.bed -m 3 -f 180 -c 1 > 180.log 2> 180.std
#python /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/TEMP/TEMP2GFF.py --insertion /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/HEG4_2.3.MSU7_BWA.insertion.refined.bp.summary 

#bash /rhome/cjinfeng/software/tools/TEMP/scripts/TEMP_Insertion_qsub.sh -i /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/bam/FC70.1.MSU7_BWA.bam -s /rhome/cjinfeng/software/tools/TEMP/scripts -r /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping_pogo.fa -t /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa.RepeatMasker.out.bed -m 3 -f 400 -c 1 > 400.log 2> 400.std
#python /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/TEMP/TEMP2GFF.py --insertion /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/HEG4_FC70_1.MSU7_BWA.insertion.refined.bp.summary

#python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE.py --te_fasta /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping.fa --genome_fasta /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa --bam /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/bam/HEG4.2.4.MSU7_BWA.bam --outdir /rhome/cjinfeng/BigData/00.RD/RILs/HEG4_mPing/RelocaTEi_HEG4.2.4 --reference_ins /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa.RepeatMasker.out --size 180


echo "merge HEG4 mPing"
bedtools intersect -wa -a RelocaTEi_HEG4.2.3/repeat/results/ALL.all_nonref_insert.characTErized.gff -b HEG4.mping.non-ref.gff > HEG4.mping.non-ref.strand.gff
bedtools window -w 5 -v -a RelocaTEi_HEG4.2.3/repeat/results/ALL.all_nonref_insert.characTErized.gff -b HEG4.mping.non-ref.gff | grep "homo" > HEG4.2.3.only.mping.temp.gff
bedtools window -w 5 -v -a HEG4.2.3.only.mping.temp.gff -b HEG4.Ping.non-ref.gff > HEG4.2.3.only.mping.gff
rm HEG4.2.3.only.mping.temp.gff
bedtools window -w 5 -v -a RelocaTEi_HEG4_P/repeat/results/ALL.all_nonref_insert.characTErized.gff -b HEG4.mping.non-ref.gff | grep "homo" > HEG4.P.only.mping.temp.gff
bedtools window -w 5 -v -a HEG4.P.only.mping.temp.gff -b HEG4.Ping.non-ref.gff > HEG4.P.only.mping.gff
rm HEG4.P.only.mping.temp.gff
bedtools intersect -v -a HEG4.P.only.mping.gff -b HEG4.2.3.only.mping.gff > HEG4.P.only.mping.toadd.gff
##include 7 ping insertion in HEG4
#cat HEG4.mping.non-ref.gff HEG4.P.only.mping.toadd.gff HEG4.2.3.only.mping.gff HEG4.Ping.non-ref.gff > HEG4.ALL.mping.non-ref.gff
cat HEG4.mping.non-ref.strand.gff HEG4.P.only.mping.toadd.gff HEG4.2.3.only.mping.gff HEG4.Ping.non-ref.gff > HEG4.ALL.mping.non-ref.gff

grep "Shared" HEG4.mping.all_inserts.gff > HEG4.mping.shared.gff
grep "only" HEG4.mping.all_inserts.gff > HEG4.mping.ref_only.gff

echo "Two library mPing: HEG4.2.3 and HEG4.P"
cat RelocaTEi_HEG4.2.3/repeat/results/ALL.all_nonref_insert.characTErized.gff RelocaTEi_HEG4_P/repeat/results/ALL.all_nonref_insert.characTErized.gff > HEG4.2libary_mPing.characTErized.gff
bedtools window -w 10 -a HEG4.ALL.mping.non-ref.AF0.1.gff -b HEG4.2libary_mPing.characTErized.gff -u | wc -l

echo "Landrace mPing"
mkdir Landrace_PNAS
cp -R ~/Rice/Rice_population_sequence/Rice_3000/Manuscript/figure/Figure_Landrace_call/*.gff Landrace_PNAS/
cp -R ~/Rice/Rice_population_sequence/Rice_3000/Manuscript/figure/Figure_Landrace_call/Landrace.call.txt Landrace_PNAS/

echo "use HEG4.P only"
#call mPing in HEG4.P
ln -s ~/Rice/RIL/genotypes_Parents/HEG4_P/HEG4_P_ATCACG_FC193L6.recal.bam HEG4_P.bam
ln -s ~/Rice/RIL/genotypes_Parents/HEG4_P/HEG4_P_ATCACG_FC193L6.recal.bai HEG4_P.bam.bai
python ReNameSRA_RelocaTEi_mPing_RelocaTE2.py --input RelocaTE2_fq --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa > log 2>&1 &
#add --bam and rerun
perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 5 --lines 1 --interval 120 --task 12 --mem 40G --time 100:00:00 --convert no RelocaTE2_fq_RelocaTEi.run.sh &
#parental mPing in HEG4.P
ln -s ../Prepare0_mPing_frequency/RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mPing.frequency ./
#python mPing_allele_frequency.py --frequency RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mPing.frequency --gff RelocaTEi_HEG4_P/repeat/results/ALL.all_nonref_insert.characTErized.gff > HEG4.P.mping.non-ref.allele.frq
#awk '$2!~/NA/' HEG4.P.mping.non-ref.allele.frq| awk '{print $1"\t"$2"\t"$2/272}' | sort -k2,2n > HEG4.P.mping.non-ref.allele.frq.sorted
python mPing_allele_frequency.py --frequency RILs_ALL_fastq_correct_merged_duplicate_RelocaTEi.CombinedGFF.characterized.clean.mPing.frequency --gff RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/ALL.all_nonref_insert.characTErized.gff > HEG4.P.RelocaTE2.mping.non-ref.allele.frq
awk '$2!~/NA/' HEG4.P.RelocaTE2.mping.non-ref.allele.frq | awk '{print $1"\t"$2"\t"$2/272}' | sort -k2,2n > HEG4.P.RelocaTE2.mping.non-ref.allele.frq.sorted
#difference between RelocaTE2 and previous version
#12
bedtools window -w 10 -a RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/ALL.all_nonref_insert.characTErized.gff -b RelocaTEi_HEG4_P/repeat/results/ALL.all_nonref_insert.characTErized.gff -v | wc -l
#4
bedtools window -w 10 -b RelocaTE2_fq_RelocaTEi/HEG4P_RelocaTEi/repeat/results/ALL.all_nonref_insert.characTErized.gff -a RelocaTEi_HEG4_P/repeat/results/ALL.all_nonref_insert.characTErized.gff -v | wc -l

