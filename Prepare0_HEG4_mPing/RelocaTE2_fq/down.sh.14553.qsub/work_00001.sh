#!/bin/bash
/opt/linux/centos/7.x/x86_64/pkgs/aspera/3.3.3/bin/ascp -i /opt/linux/centos/7.x/x86_64/pkgs/aspera/3.3.3/etc/asperaweb_id_dsa.openssh -k 1 -T -l20m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR631/SRR631735/SRR631735.sra /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Manuscript_Preparation/Prepare0_HEG4_mPing/RelocaTE2_fq/Citrus_RNAseq/SRR631735; echo This-Work-is-Completed!