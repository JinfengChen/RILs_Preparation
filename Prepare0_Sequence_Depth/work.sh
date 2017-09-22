echo "depth"
sed 's/GN//' RILs_ALL_bam_correct_merged.summary | sort -k1n > RILs_ALL_bam_correct_merged.sorted.summary

echo "coverage"
cp ~/BigData/00.RD/RILs/QTL_pipe/bin/RILs_ALL_272line/NB.RILs.dbSNP.SNPs.RILs ./
python NeedCare.py > RILs_ALL_272.needcare.txt &

echo "Generate Depth.table"
python SeqDepth.py --input RILs_ALL_bam_correct_merged.sorted.summary

echo "Draw Depth curve"
cat Depth.R | R --slave

echo "Draw Depth vs Genotype% plot"
cat DepthvsGT.R | R --slave

