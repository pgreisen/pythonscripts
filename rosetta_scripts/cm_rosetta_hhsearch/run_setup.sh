for i in $(seq 1 100);
do
    if [ -d "$i" ]; then
	cd $i;
	fasta=`find *fasta`;
	hhsearch=`find *_pdb.hhr`;
	python ~/bin/get_top_n_sequence.py $hhsearch 10;
	python2.7 ~/bin/setup_RosettaCM.py --fasta $fasta --alignment reduced_hhsearch.hhr --alignment_format hhsearch --rosetta_bin ~/bin/rosetta_bin --build static > /dev/null & echo "setup dir $i";
	cd ..;
    fi
done
