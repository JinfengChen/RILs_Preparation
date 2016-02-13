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

