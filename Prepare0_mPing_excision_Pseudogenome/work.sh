#get flanking 20bp sequence to valid pseudogenome
bedtools flank -i Parent.ALL.mPing.Ref_Shared.gff -g MSU_r7.fa.fai -b 20 > Parent.ALL.mPing.Ref_Shared.flank20bp.gff
bedtools getfasta -fi MSU_r7.fa -bed Parent.ALL.mPing.Ref_Shared.flank20bp.gff -fo Parent.ALL.mPing.Ref_Shared.flank20bp.fa
#generate pseudogenome and valid the 20bp sequence by print
python Pseudo_TEinsertion_Genome_Ref.py --gff Parent.ALL.mPing.Ref_Shared.gff --genome MSU_r7.fa --project Parent.Pseudo_mPing.Ref_Shared
#valid by extracing 20 bp sequence from gff
bedtools flank -i Parent.Pseudo_mPing.Ref_Shared.gff -g Parent.Pseudo_mPing.Ref_Shared.fa.fai -b 20 > Parent.Pseudo_mPing.Ref_Shared.flank20bp.gff
bedtools getfasta -fi Parent.Pseudo_mPing.Ref_Shared.fa -bed Parent.Pseudo_mPing.Ref_Shared.gff -fo Parent.Pseudo_mPing.Ref_Shared.TSD.fa
bedtools getfasta -fi Parent.Pseudo_mPing.Ref_Shared.fa -bed Parent.Pseudo_mPing.Ref_Shared.flank20bp.gff -fo Parent.Pseudo_mPing.Ref_Shared.flank20bp.fa
#check on with 1 bp less: Chr10	HEG4	transposable_element_insertion_site	21716391	21716819 
awk '{print $5-$4,$_}' Parent.ALL.mPing.Ref_Shared.gff | less -S
bedtools getfasta -fi MSU_r7.fa -bed Parent.ALL.mPing.Ref_Shared.gff -fo Parent.ALL.mPing.Ref_Shared.fa

