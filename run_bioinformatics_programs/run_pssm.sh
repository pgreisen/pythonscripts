path=/home/ubuntu/ncbi-blast-2.7.1+/bin
database=/home/ubuntu/uniref90
$path/psiblast -query $1 -db $database/uniref90_w_index.fasta -out_pssm my_protein.ckp -evalue 0.01 -out_ascii_pssm ascii_mtx_file -out output_file -num_iterations 3 