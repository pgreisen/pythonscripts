rm 20160428_human.fasta;
for i in *V_HUMAN_1.fasta; 
do 
    perl /z/insulin/users/pjug/Programs/igblast/edit_imgt_file.pl $i >> 20160428_human_VH.fasta; 
done
/z/insulin/users/pjug/Programs/igblast/bin/makeblastdb -in 20160428_human_VH.fasta -out 20160428_human_VH.db -dbtype prot -logfile 20160428_human.log -parse_seqids

for i in *J_HUMAN_1.fasta; 
do 
    perl /z/insulin/users/pjug/Programs/igblast/edit_imgt_file.pl $i >> 20160428_human_VJ.fasta; 
done
/z/insulin/users/pjug/Programs/igblast/bin/makeblastdb -in 20160428_human_VJ.fasta -out 20160428_human_VJ.db -dbtype prot -logfile 20160428_human.log -parse_seqids

for i in *D_HUMAN_1.fasta; 
do 
    perl /z/insulin/users/pjug/Programs/igblast/edit_imgt_file.pl $i >> 20160428_human_VD.fasta; 
done
/z/insulin/users/pjug/Programs/igblast/bin/makeblastdb -in 20160428_human_VD.fasta -out 20160428_human_VD.db -dbtype prot -logfile 20160428_human.log -parse_seqids
