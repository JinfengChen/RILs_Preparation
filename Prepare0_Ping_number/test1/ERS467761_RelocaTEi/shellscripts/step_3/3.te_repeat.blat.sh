/opt/linux/centos/7.x/x86_64/pkgs/blat/35/bin/blat -minScore=10 -tileSize=7 /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/ping.fa /bigdata/wesslerlab/shared/Rice/Rice_population_sequence/Rice_3000/CAAS/Japonica_fastq/ERS467761/ERR615110_2.fa /bigdata/wesslerlab/shared/Rice/Rice_population_sequence/Rice_3000/CAAS/Japonica_fastq_RelocaTEi/ERS467761_RelocaTEi/repeat/blat_output/ERR615110_2.te_repeat.blatout 1>>/bigdata/wesslerlab/shared/Rice/Rice_population_sequence/Rice_3000/CAAS/Japonica_fastq_RelocaTEi/ERS467761_RelocaTEi/repeat/blat_output/blat.out 2>>/bigdata/wesslerlab/shared/Rice/Rice_population_sequence/Rice_3000/CAAS/Japonica_fastq_RelocaTEi/ERS467761_RelocaTEi/repeat/blat_output/blat.out
python /rhome/cjinfeng/BigData/00.RD/RelocaTE2_mPing/scripts/relocaTE_trim.py /bigdata/wesslerlab/shared/Rice/Rice_population_sequence/Rice_3000/CAAS/Japonica_fastq_RelocaTEi/ERS467761_RelocaTEi/repeat/blat_output/ERR615110_2.te_repeat.blatout /bigdata/wesslerlab/shared/Rice/Rice_population_sequence/Rice_3000/CAAS/Japonica_fastq/ERS467761/ERR615110_2.fastq.gz 10 1 > /bigdata/wesslerlab/shared/Rice/Rice_population_sequence/Rice_3000/CAAS/Japonica_fastq_RelocaTEi/ERS467761_RelocaTEi/repeat/flanking_seq/ERR615110_2.te_repeat.flankingReads.fq
