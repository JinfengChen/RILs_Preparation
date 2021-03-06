#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

#Nuc ch1
#python TSSprofile_DHS_R.py ../input/DHS.Chr1.unique.bam ../input/MSU7.gene.exon_number.Chr1.gtf
#python TTSprofile_DHS_R.py ../input/DHS.Chr1.unique.bam ../input/MSU7.gene.exon_number.Chr1.gtf
#Nuc genome
python TSSprofile_DHS_R.py ../input/DHS.unique.bam ../input/MSU7.gene.exon_number.gtf
python TTSprofile_DHS_R.py ../input/DHS.unique.bam ../input/MSU7.gene.exon_number.gtf

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

