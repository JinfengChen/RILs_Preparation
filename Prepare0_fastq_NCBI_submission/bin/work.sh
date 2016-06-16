#NCBI biosample templetate
#../input/Plant.1.0.tsv
#NCBI SRA meta
#../input/SRA_metadata.tsv
echo "Generate batch biosample for 272 RILs"
python NCBI_Biosample_batch.py --input RILs_ALL_fastq_correct_merged.fastq.list
echo "Generate batch biosample for 272 RILs"
python NCBI_SRA_batch.py --input RILs_ALL_fastq_correct_merged.fastq.list

echo "check file consistence"
ls ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_fastq_correct_merged_duplicate/RIL*/*_p1.fq.gz > dupli.list
ls ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_fastq_correct_merged/RIL*/*_p1.fq.gz > merged.list
grep "^RIL" RILs_272.NCBI_SRA.tsv | cut -f13> sra.list
python listdiff.py merged.list sra.list


echo "update to NCBI"
ssh m05
mkdir /scratch/RILs_submission 
cd /scratch/RILs_submission
ln -s ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_fastq_correct_merged/RIL*/*.gz ./
bash md5.sh > md5.log 2>&1
ls *.gz | awk '{print "/opt/aspera/3.3.3/bin/ascp -i /opt/aspera/3.3.3/etc/asperaweb_id_dsa.openssh -QT -l100m -k1 -d " $1 " subasp@upload.ncbi.nlm.nih.gov:uploads/jinfeng.chen@ucr.edu_yWAL07Qi/PRJNA316308/"}' > upload.sh
bash upload.sh > upload.log 2>&1
bash upload_md5sum.sh > upload_md5sum.log 2>&1 &

