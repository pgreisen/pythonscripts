rm *fasta;
for i in 1;
do 
    python ~/pythonscripts/develop/generate_variants/combine_positions.py muts.txt $i cbdas.pep 
done
