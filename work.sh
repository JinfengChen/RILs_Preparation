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
