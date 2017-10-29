cd ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Pseudogenome
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/HEG4.2.3.RelocaTE2.mping.non-ref.AF0.1.gff ./
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Identification/input/Parent.ALL.mPing.Ref_Shared.gff ./
python Pseudo_TEinsertion_Genome.py --repeat mPing_Ping_Pong.fa --gff HEG4.2.3.RelocaTE2.mping.non-ref.AF0.1.gff --genome MSU_r7.fa --project MSU_r7.Pseudo_mPing_415

cd /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Identification/input
ln -s /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Pseudogenome/MSU_r7.Pseudo_mPing_415.gff ./
ln -s /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Pseudogenome/HEG4.2.3.RelocaTE2.mping.non-ref.AF0.1.gff Parent.ALL.mPing.415.gff
grep -e "Shared" -e "Reference" -e "TE=ping" Parent.ALL.mPing.gff > Parent.ALL.mPing.others.gff
grep -v "TE=ping" Parent.ALL.mPing.others.gff > Parent.ALL.mPing.Ref_Shared.gff
cat ../../Prepare0_mPing_excision_Pseudogenome/HEG4.2.3.RelocaTE2.mping.non-ref.AF0.1.gff Parent.ALL.mPing.others.gff | sort -k1,1 -k4,4n > Parent.ALL.mPing.415_plus_other.gff
python PrepareFastq.py --bam ../input/RILs_ALL_bam_correct_merged > log 2>&1 &
python runRIL_bwa.py --input ../input/RILs_ALL_unmapped_mping_fastq > log 2>&1 &

#fix RIL105 duplicate read bug
cd ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_excision_Identification/input/RILs_ALL_unmapped_mping_fastq/RIL105
#/rhome/cjinfeng/BigData/software/FastUniq/source/fastuniq -i file.list -o RIL105_uniq_1.fq -p RIL105_uniq_2.fq
sbatch uniq.sh
mv RIL105_uniq_1.fq RIL105_1.fq
mv RIL105_uniq_2.fq RIL105_2.fq

cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_272line/MPR.geno.data ../input/
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_272line/MPR.geno.bin ../input/
mv RIL*.bam ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing/
mv RIL*.bai ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing/
mv RIL*.dupli ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing/

#python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam_correct_merged --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing --gff_ref ../input/Parent.ALL.mPing.415.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing_415.gff > log 2>&1 &
sbatch run_excision.sh
python Ping_number_RILs.High_exicison.py --csv mPing_boundary_415_nonref_mPing --ping_code ../../Prepare0_mPing_calls/RIL272_RelocaTEi.Jinfeng_Lulu.ping_code.table.txt
mv mPing_boundary_mPing_GT_Ping_code mPing_boundary_415_nonref_mPing_GT_Ping_code
python Sum_excision_distance.py --dir mPing_boundary_415_nonref_mPing_GT_Ping_code --distance ../../Prepare0_mPing_distance/mPing_dist2.50Mb.list.sorted --blacklist Bam.Core.blacklist --project mPing_boundary_415_nonref_mPing.linked_50Mb_debug2 > log 2>&1 &

