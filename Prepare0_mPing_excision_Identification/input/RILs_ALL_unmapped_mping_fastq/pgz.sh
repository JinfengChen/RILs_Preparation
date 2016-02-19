#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -l mem=20gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR


for i in `ls RIL*/*_1.fq | sed 's/@//'`
do
   echo $i
   if [ ! -e $i.gz ]; then
      pigz $i -p 16
      chmod 664 $i.gz
      #gzip $i
   fi
done

echo "Done"
