
for i in `ls *.gz | sed 's/@//'`
do
   echo $i
   md5sum $i >> md5sum.list
done


