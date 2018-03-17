grep "Shared" HEG4_2.3.mPing.20X.mping.all_inserts.gff > HEG4.Shared_mPing.gff
grep "Ref" HEG4_2.3.mPing.20X.mping.all_inserts.gff > HEG4.RefOnly_mPing.gff
cat HEG4.ALL.mping.non-ref.AF0.1.gff HEG4.RefOnly_mPing.gff HEG4.Shared_mPing.gff > HEG4_NB.mPing.gff
cat HEG4.RefOnly_mPing.gff HEG4.Shared_mPing.gff > HEG4.NB_mPing.gff


echo "Ref_Only: use distance in NB as they present only in NB"
perl mPing_dist_ref.pl --input HEG4.NB_mPing.gff | grep "ref_only" | sort -k1,1 -k2,2n > Distance.Refonly.txt
echo "Shared: use distance either in HEG4 or Nipponbare whoever is short, because they present in both"
perl mPing_dist_ref.pl --input HEG4_NB.mPing.gff | grep "share" | sort -k1,1n -k2,2n > Distance.Shared.txt


echo "HEG4 and NB 466 mPing"
ln -s ~/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_mPing_calls/HEG4_NB_ALL.mPing.gff ./
perl mPing_dist.pl --input HEG4_NB_ALL.mPing.gff
sort -k3,3n -k1,1n -k2,2n mPing_dist2.50Mb.list | uniq > mPing_dist2.50Mb.list.sorted

