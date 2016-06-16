echo "matrix of genotype high excision mPing pairs for all the RILs"
#python mPing_unknown_gt_pairs_RILs.py --input High_excision_mPing.pairs --matrix ../input/High_excision_csv_Ping
python mPing_unknown_gt_pairs_RILs.py --input High_excision_mPing.pairs --matrix ../input/mPing_boundary_mPing_RILs_GT_Ping_code

echo "RILs for PCR, which are not NB, not in blacklist of problem RILs"
#need to manual edit in excel during the process and also label these we already have test for high excision mping
python RILs_need.py > High_excision_mPing.needDNA

echo "get bam file to check"
python get_excision_bam.py --input High_excision_mPing.draw.SV --output High_excision_mPing_draw_SV_bam
tar -zcvf High_excision_mPing_draw_SV_bam.tar.gz High_excision_mPing_draw_SV_bam/**
